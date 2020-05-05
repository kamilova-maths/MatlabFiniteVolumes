% Clear previous files
clear all 
clc

% Parameters shared with other routines (I'm not very sure about this,but
% let's give it a try)
global Pe Bi tha N K gamma P0 St T L D uf 

% Define parameters
% Data to use 

rho= 1.8*10^3; %Bergstrom ; 
g = 10; 
c= 900; % Fitt and Howell
L=7; % Temperature Profiles in Soderberg Electrodes
uc = 10^-5; %Bergstrom approximation
R0=0.5;
R1=1; 
k = 3; 
Qc=15000; %Taken very vaguely from Temperature profiles in Soderberg electrodes. 
mu0=10^10; % given by Bjornar at a reference temperature
T_a=343; 
h = 7; 

%Defining non-dimensional parameters
% Peclet number
Pe = (rho*c*uc*L)/(k);
epsilon=R1/L;
St=(rho*g*L^2)/(uc*mu0);
P0 = (10000*L)/((R1^2)*uc*mu0);

DeltaT = (Qc*L)/(rho*c*uc);
Bi= ((L^2)*h)/(k*R1); 
tha = 0.005; 
D = (R0^2)/(R1^2); 

gamma = 40; 

%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
eps = 1e-2;

% Calculating the initial conditions as a solution of the steady state
% problem 
N=1000; K=300;
% end of the domain
T = 1; L=1 ;
dx = 1/K;

% We add the heaviside with H=1, and we remove it with H=0. 
H=1;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;

% We try with tha=0
[P, A0, J0, th0, ~] = InitialConditionsSteady(3*K,gamma,Q,x1,x2,eps,St,tha,Bi,Pe,P0,R0,L,H,plt);

% Find lam0, then resize both sides with an interpolation (only necessary
% to do this for the steady state, the conditions on the rest are much
% nicer because of our rescalings 
I =    find(A0>0.999,1,'first');
x = linspace(0,1,3*K);
lam0 = x(I);

A0top = interp1(x(1:I),A0(1:I),linspace(0,x(I),K),'spline'); 
th0top = interp1(x(1:I),th0(1:I),linspace(0,x(I),K),'spline');

%A0bot = ones(size(A0top));
th0bot = interp1(x(I:end),th0(I:end),linspace(x(I),1,K),'spline');


%Initial conditions given by A0, th0
A0=A0top';
th0=th0top'; 

% Obtain initial condition for u
% viscosity
%mu = @(th) exp(-gamma*th); % Perhaps discretise this as well, but put a pin in that
% solve for u
temp1 = [0;0;th0];
temp = (temp1(1:end-1)+temp1(2:end) ) /2; 
tiph= 3*exp(-gamma*temp); % evaluation  
%tiph = 1.6*ones(size(temp));
tmpA = [ D; A0  ];           % add ghost node to A

% I = find(A0>0.999,1,'first');
% x=0:dx:L;

% Note that the steady state we have used has gamma = 40, which is
% different to what we are doing for constant theta ... This is why it is
% important that I write the following code leaving the necessary blanks
% for the mu(theta) portion of the code 

uf= 1/tmpA(end);
Dx2u = spdiags( [ tmpA(2:end).*tiph(2:end), -(tmpA(1:end-1).*tiph(1:end-1)+tmpA(2:end).*tiph(2:end)), tmpA(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], K, K );
Dx2u(1,2) = Dx2u(1,2) + tmpA(1).*tiph(1) / dx^2;              % include effect from Neumann BC
fu   = - (lam0^2)*St*( tmpA(1:end-1) + tmpA(2:end) )/ 2;
fu(1) = fu(1) -2*tiph(1)*P0(1)/(3*dx); 
fu(end)= fu(end) - uf   * tmpA(end)* tiph(end)/dx^2;
u0 = Dx2u\fu;

u0top = u0; 
%u0bot = ones(size(u0top));


y0(1:K) = A0top;
y0(1+K:2*K) = u0top;
y0(2*K+1:3*K) = th0top; 

y0(3*K+1:4*K) = lam0.*ones(size(u0));

% y0(4*K+1:5*K) = A0bot;
% y0(5*K+1:6*K) = u0bot;
y0(4*K+1:5*K) = th0bot;


% Independent variable for ODE integration 
tout = linspace(0,T,N);

% ODE integration 
reltol = 1.0e-04; abstol = 1.0e-04;
options = odeset('RelTol',reltol,'AbsTol',abstol);
%[t,y] = ode15s(@coupledPdeNoTemp,tout,y0);
[t,y] = ode15s(@coupledPde,tout,y0); 

Atop  = y(:,1:K);
utop  = y(:,1+K:2*K);
th = [y(:,2*K+1:3*K),y(:,4*K+1:5*K)];

lam   = y(:,3*K+1:4*K);
% 
% Abot  = y(:,4*K+1:5*K); % This should just be zeros


% We try again but with the achieved steady state from the previous method 

%% FOR NOW, IGNORE HOW WE OBTAIN THE SOLUTION TO U 
%% Solve for u at next time step (laplacian)
%     temp = 3*mu([0;0;th(:,i)]);     % add ghost node to th
%     tiph = ( temp(1:end-1) + temp(2:end) ) / 2;
Aforu = Atop';
thforu = y(:,2*K+1:3*K)'; 
u = zeros(size(Aforu));
u(:,1) = u0; 
for i=2:N
   
    %% Solve for u at next time step (laplacian)
%     temp = 3*mu([0;0;th(:,i)]);     % add ghost node to th
%     tiph = ( temp(1:end-1) + temp(2:end) ) / 2;
    
    % Technically can do this in one line, but I can't figure out how -
    % this gives same results as above so we don't bother changing it
    temp1 = [0;0;thforu(:,i)];
    temp = (temp1(1:end-1)+temp1(2:end) ) /2; 
    tiph= 3*exp(-gamma*temp); 
    
    Ag  = [ D; Aforu(:,i)  ];           % add ghost node to A

    Dx2u = spdiags( [ Ag(2:end).*tiph(2:end), -(Ag(1:end-1).*tiph(1:end-1)+Ag(2:end).*tiph(2:end)), Ag(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], K, K );
    Dx2u(1,2) = Dx2u(1,2) + Ag(1).*tiph(1) / dx^2;              % include effect from Neumann BC
    fu   = - St*( Ag(1:end-1) + Ag(2:end) )/ 2;

    if length(P0)==1
    fu(1) = fu(1) -2*tiph(1)*P0/(3*dx);  % include derivative (again, Neumann BC)
    else
    fu(1) = fu(1) -2*tiph(1)*P0(i)/(3*dx);  % include derivative (again, Neumann BC)
    end
    
    fu(end)= fu(end) - uf   * Ag(end)* tiph(end)/dx^2;
    u(:,i) = Dx2u\fu;
 
    
end

%thfullT = [ zeros(N,1), ...
    %   thtop          ];
   
thfull = [ zeros(N,1), ...
       th         ];
  
Afull  = [ D*ones(N,1), ...
       Atop           ];
   
ufull  = [ u' , ...
       uf.*ones(N,1) ];
   
dx = 1/K;

x = (0:dx:L)';


figure;
surf(t,x,Afull','LineStyle','none')
xlabel('$t$')
ylabel('$X$')
title('$A$ with MOL')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

figure;
surf(t,linspace(0,1,2*K+1),thfull','LineStyle','none')
xlabel('$t$')
ylabel('$x$')
title('$\theta$ with MOL')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

figure;
surf(t,x,ufull','LineStyle','none')
xlabel('$t$')
ylabel('$X$')
title('$u$ with MOL')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

plot(t,lam(:,end))
xlabel('$t$')
ylabel('$\lambda(t)$')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

