% Clear previous files
clear all 
clc

% Parameters shared with other routines (I'm not very sure about this,but
% let's give it a try)
global Pe Bi tha N K gamma P0 St T L D uf x1 x2 Q

% Define parameters
% Data to use 

rho= 1.8*10^3; %Bergstrom ; 
g = 10; 
c= 900; % Fitt and Howell
Ld=7; % Temperature Profiles in Soderberg Electrodes
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
Pe = (rho*c*uc*Ld)/(k);
epsilon=R1/Ld;
St=(rho*g*Ld^2)/(uc*mu0);
P0 = (10000*Ld)/((R1^2)*uc*mu0);

DeltaT = (Qc*Ld)/(rho*c*uc);
Bi= ((Ld^2)*h)/(k*R1); 
tha = 0.005; 
D = (R0^2)/(R1^2); 

gamma = 20; 

%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
eps = 1e-2;

% Calculating the initial conditions as a solution of the steady state
% problem 
N=1000; K=600;
% end of the domain
T = 1; L=1.5 ;
dx = L/K;
uf = 1; 

% We add the heaviside with H=1, and we remove it with H=0. 
H=1;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;

% We try with tha=0
[P, A0, J0, th0, ~] = InitialConditionsSteady(4*K,eps,H,plt);
%pause
%close all 


% Find lam0, then resize both sides with an interpolation (only necessary
% to do this for the steady state, the conditions on the rest are much
% nicer because of our rescalings 
I =    find(A0>0.99999,1,'first');
x = linspace(0,L,4*K);
lam0 = x(I);

A0top = interp1(x(1:I),A0(1:I),linspace(0,x(I),K),'spline'); 
th0top = interp1(x(1:I),th0(1:I),linspace(0,x(I),K),'spline');

%A0bot = ones(size(A0top));
th0bot = interp1(x(I:end),th0(I:end),linspace(x(I),L,K),'spline');


%Initial conditions given by A0, th0
A0=A0top';
th0=th0top'; 
th0bot = flip(th0bot)'; 

% Obtain initial condition for u
% viscosity
%mu = @(th) exp(-gamma*th); % Perhaps discretise this as well, but put a pin in that
% solve for u

% I = find(A0>0.999,1,'first');
% x=0:dx:L;

% Note that the steady state we have used has gamma = 40, which is
% different to what we are doing for constant theta ... This is why it is
% important that I write the following code leaving the necessary blanks
% for the mu(theta) portion of the code 


%u0bot = ones(size(u0top));
y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = lam0.*ones(size(A0));  

y0(3*K+1:4*K) = th0bot;

% Independent variable for ODE integration 
tout = linspace(0,T,N);

% ODE integration 
reltol = 1.0e-04; abstol = 1.0e-04;
options = odeset('RelTol',reltol,'AbsTol',abstol);
%[t,y] = ode15s(@coupledPdeNoTemp,tout,y0);
tic
[t,y] = ode15s(@coupledPde,tout,y0); 
% [t,y] = ode15s(@coupledPde2,tout,y0); 
toc
%pause
A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

lam   = y(:,2*K+1:3*K);

phi = y(:,3*K+1:4*K); 

u = zeros(size(A));
u0 = usolutionNoFB(A0,th0);  % Do I even need this guy? 
u(1,:) = u0; 
for i=2:N
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i,end),L);   
    
end
disp('Completed Round 1')

thtop = [ zeros(N,1), ...
       th         ];
  
Atop  = [ D*ones(N,1), ...
       A           ];
   
tmpA =  (Atop + [Atop(:,2:end), ones(N,1)])/2; % extract A at the edges 
Abot = ones(N,K); 
Afull =    [tmpA, Abot];  
%ufull  = [ u , ...
 %      uf.*ones(N,1) ];
   
dx = 1/K;

x = (0:dx:1)';
[X, T1] = meshgrid(x,t);
Xresc1 = [lam(:,1), lam].*X; 
[X, T2] = meshgrid((dx:dx:1)',t);
Xresc2 = -X.*(L-lam)+L;

numel=10;
datamat = [[Xresc1(1,:)'; flip(Xresc2(1,:))'], [thtop(1,:)'; flip(phi(1,:))']];

figure; 
for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),thtop(i,:))
    hold on
    plot(Xresc2(i,:),phi(i,:)) 
   % pause
    datamat = [datamat, [[Xresc1(i,:)';flip(Xresc2(i,:))'], [thtop(i,:)'; flip(phi(i,:))']]];
end
set(gca,'TickLabelInterpreter','latex','fontsize',13)

csvwrite('ThetaDiscreteTimesteps.csv',datamat); 
return

% Round 2 
A0  = A(end,:);
%u0 = u(end,:);
th0 = th(end,:);
th0bot = flip( phi(end,:) ); 

y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = lam(end,:);  

y0(3*K+1:4*K) = th0bot;
% 

tic
[t,y] = ode15s(@coupledPde,tout,y0); 
toc
A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

lam   = y(:,2*K+1:3*K);

phi = y(:,3*K+1:4*K); 

u = zeros(size(A));
u0 = usolutionNoFB(A0',th0');  % Do I even need this guy? 
u(1,:) = u0; 
for i=2:N
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i,:)', L);   
    
end


% Abot  = y(:,4*K+1:5*K); % This should just be zeros


% We try again but with the achieved steady state from the previous method 

%% FOR NOW, IGNORE HOW WE OBTAIN THE SOLUTION TO U 
%% Solve for u at next time step (laplacian)
%     temp = 3*mu([0;0;th(:,i)]);     % add ghost node to th
%     tiph = ( temp(1:end-1) + temp(2:end) ) / 2;

%thfullT = [ zeros(N,1), ...
    %   thtop          ];
   
thtop = [ zeros(N,1), ...
       th         ];
  
Atop  = [ D*ones(N,1), ...
       A           ];
   
tmpA =  (Atop + [Atop(:,2:end), ones(N,1)])/2; % extract A at the edges 
Abot = ones(N,K); 
Afull =    [tmpA, Abot];  
ufull  = [ u , ...
       uf.*ones(N,1) ];
   
dx = L/K;

x = (0:dx:L)';
[X, T1] = meshgrid(x,t);
Xresc1 = [lam(:,1), lam].*X; 
[X, T2] = meshgrid((dx:dx:L)',t);
Xresc2 = -X.*(L-lam)+L;
set(gca,'TickLabelInterpreter','latex','fontsize',13)

% figure;
% surf(T1,Xresc1,thtop,'LineStyle','none')
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)


% figure; 
% surf(T2,Xresc2,phi,'LineStyle','none')
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% 
% figure; 
% surf(T1,Xresc1,thtop,'LineStyle','none')
% hold on 
% surf(T2,Xresc2,phi,'LineStyle','none')
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)

numel=20;
datamat = [[Xresc1(1,:)'; flip(Xresc2(1,:))'], [thtop(1,:)'; flip(phi(1,:))']];

figure; 
for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),thtop(i,:))
    hold on
    plot(Xresc2(i,:),phi(i,:)) 
    pause
    datamat = [datamat, [[Xresc1(i,:)';flip(Xresc2(i,:))'], [thtop(i,:)'; flip(phi(i,:))']]];
end
%csvwrite('ThetaDiscreteTimesteps.csv',datamat); 



numel=10;
datamat = [[Xresc1(1,:)'; flip(Xresc2(1,:))'], [tmpA(1,:)'; flip(Abot(1,:))']];
figure; 
for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),tmpA(i,:))
    hold on
    plot(Xresc2(i,:),Abot(i,:))
    datamat = [datamat, [[Xresc1(i,:)';flip(Xresc2(i,:))'], [tmpA(i,:)'; flip(Abot(i,:))']]];
end
% 
% %csvwrite('ADiscreteTimesteps.csv',datamat); 
% 
% numel=10;
% datamat = [[Xresc1(1,:)'; flip(Xresc2(1,:))'], [ufull(1,:)'; flip(Abot(1,:))']];
% figure; 
% for i = N/numel:(N/numel):N
%     plot(Xresc1(i,:),ufull(i,:))
%     hold on
%     plot(Xresc2(i,:),Abot(i,:))
%     datamat = [datamat, [[Xresc1(i,:)';flip(Xresc2(i,:))'], [ufull(i,:)'; flip(Abot(i,:))']]];
% end
% csvwrite('uDiscreteTimesteps.csv',datamat);
% 
% this is a file that has two columns per timestep, where the first column
% is the x values and the second one the theta values. In the x and theta
% values, we have appended both sides onto one vector, so Xcol = Xresc1,
% Xresc2, and similarly for theta 



return 
figure;
% surf(t,x.*[lam(1,1); lam(:,end)],Atop','LineStyle','none')
% xlabel('$t$')
% ylabel('$X$')
% title('$A$ with MOL','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% figure;
% surf(t,x,thfull','LineStyle','none')
% xlabel('$t$')
% ylabel('$x$')
% title('$\theta$ with MOL','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)

figure;
surf(t,x,ufull','LineStyle','none')
xlabel('$t$')
ylabel('$X$')
title('$u$ with MOL','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

plot(t,lam(:,end))
xlabel('$t$','Interpreter','latex')
ylabel('$\lambda(t)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

datamat = [[t(1,end); t(10:10:end)],[lam(1,end); lam(10:10:end,end)]];

csvwrite('lambdavst.csv',datamat);

figure;
%surf(T1,Xresc1,Atop,'LineStyle','none')
surf(T3,Xresc3,Afull,'LineStyle','none')

% hold on
% Abot = ones(N,K); 
%surf(T2,Xresc2,Abot,'LineStyle','none')
xlabel('$t$','Interpreter','latex')
ylabel('$x$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)



