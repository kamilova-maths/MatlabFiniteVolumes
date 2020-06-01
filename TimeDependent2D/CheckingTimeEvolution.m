% Clear previous files
clear all 
close all 
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
eps = 1e-4;

% Calculating the initial conditions as a solution of the steady state
% problem 
N=800; K=300;
% end of the domain
T = 3; L=1.5 ;
dx = L/K;
uf = 1; 

% We add the heaviside with H=1, and we remove it with H=0. 
H=1;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;


[Psteady, A0steady, J0steady, th0steady, xsteady] = InitialConditionsSteady(eps,H,plt);
% A0 and th0 are actually whatever size Matlab needs them to be, as it uses
% an adaptive mesh. It is our job to interpolate this accordingly. 
%close all 



% Find lam0, then resize both sides with an interpolation (only necessary
% to do this for the steady state, the conditions on the rest are much
% nicer because of our rescalings 
A0avg  = (A0steady(1:end-1)+A0steady(2:end))/2;   
th0avg = (th0steady(1:end-1) + th0steady(2:end))/2; 
xavg   = (xsteady(1:end-1) + xsteady(2:end))/2;
I =    find(A0avg>0.99999,1,'first');

lam0steady = xavg(I);

A0topsteady = interp1(xavg(1:I),A0avg(1:I),linspace(0,xavg(I),K),'pchip'); 
th0topsteady = interp1(xavg(1:I),th0avg(1:I),linspace(0,xavg(I),K),'pchip'); % pchip and cubic should be exactly the same

%A0bot = ones(size(A0top));
th0botsteady = interp1(xavg(I:end),th0avg(I:end),linspace(xavg(I),L,K),'pchip');

%u0steady = usolution(A0topsteady,th0topsteady,lam0steady,1);   
lam0 = lam0steady;
u0steady = usolution(A0topsteady',th0topsteady',lam0steady,1);   


A0 = (1-D).*linspace(0,1,K)'+D; 
th0 = zeros(K,1); 
% th0bot = sin(((linspace(lam0,L,K)-lam0).*pi)/(2.*L))';
th0bot = zeros(K,1); 
%th0=th0top'; 
th0bot = flip(th0bot); % We resize to Kx1 and flip the bottom part for theta,
                        % to match the requirements from coupledPde
y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = th0bot;  

y0(3*K+1) = lam0;

% Independent variable for ODE integration 
tout = linspace(0,T,N);

% ODE integration 
reltol = 1.0e-04; abstol = 1.0e-04;
options = odeset('RelTol',reltol,'AbsTol',abstol);
tic
[t,y] = ode15s(@coupledPde2,tout,y0); 
toc

A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

phi   = y(:,2*K+1:3*K);

lam = y(:,3*K+1); 

% This is the solution to u at the top. 
u = zeros(size(A));
for i=1:N
    % Here A(1,:) is A0, so this gives me the u used for the
    % calculations inside coupledPde.m
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i,end),1);   
end
disp('Completed Round 1')

thtop = [ zeros(N,1), ...
       th         ];
  
Atop  = [ D*ones(N,1), ...
       A           ];
   
tmpA =  (Atop + [Atop(:,2:end), ones(N,1)])/2; % extract A at the edges 
Abot = ones(N,K); 
Afull =    [tmpA, Abot];  
ufull  = [ u , ...
      uf.*ones(N,1) ];
   
dx = 1/K;
% we define a vector for lambda(t) [so that I don't get confused with
% matrices]. Recall that lambda only depends on t, so it is constant for
% each x


x = (0:dx:1)';
[X, T1] = meshgrid(x,t);
Xresc1 = lam.*X; 
xbar = linspace(1,0,K)'; 
[X, T2] = meshgrid(xbar,t);
Xresc2 = -X.*(L-lam)+L;

% PLOTTING THETA
numel=10;
datamat1 = [[Xresc1(1,1:5:end)'; Xresc2(1,1:5:end)'], [thtop(1,1:5:end)'; flip(phi(1,1:5:end))']];

figure; 
plot(Xresc1(1,:),thtop(1,:))
hold on 
plot(Xresc2(1,:),flip(phi(1,:)))
hold on 
for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),thtop(i,:))
    hold on
    plot(Xresc2(i,:),flip(phi(i,:))) 
    %pause
    datamat1 = [datamat1, [[Xresc1(i,1:5:end)';Xresc2(i,1:5:end)'], [thtop(i,1:5:end)'; flip(phi(i,1:5:end))']]];
end
set(gca,'TickLabelInterpreter','latex','fontsize',13)
top = max(max(abs(th(1,:)-th(end,:))))/max(max(th))
bottom = max(max(abs(phi(1,:)-phi(end,:))))/max(max(phi))
csvwrite('ThetaDiscreteTimestepsgamma20.csv',datamat1); 
hold on 
plot(Xresc1(1,:),[0,th0topsteady],'--')
hold on 
plot(Xresc2(1,:),th0botsteady,'--') 


% PLOTTING A
figure; 
plot(Xresc1(1,:),tmpA(1,:))
hold on 
plot(Xresc2(1,:),Abot(1,:))
numel=10;
datamat2 = [[Xresc1(1,1:5:end)'; Xresc2(1,1:5:end)'], [tmpA(1,1:5:end)'; flip(Abot(1,1:5:end))']];

for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),tmpA(i,:))
    hold on
    plot(Xresc2(i,:),Abot(i,:))
    datamat2 = [datamat2, [[Xresc1(i,1:5:end)'; Xresc2(i,1:5:end)'], [tmpA(i,1:5:end)'; flip(Abot(i,1:5:end))']]];
end
hold on 
plot(Xresc1(1,:),[D,A0topsteady],'--') 
%plot(Xresc1(1,:),[Atop(1,1); Atop(1,:)],'--')


csvwrite('ADiscreteTimestepsgamma20.csv',datamat2); 

% PLOTTING U
figure; 
plot(Xresc1(1,:),ufull(1,:))
hold on 
plot(Xresc2(1,:),Abot(1,:))
numel=10;
datamat3 = [[Xresc1(1,1:10:end)'; Xresc2(1,1:10:end)'], [ufull(1,1:10:end)'; flip(Abot(1,1:10:end))']];

for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),ufull(i,:))
    hold on
    plot(Xresc2(i,:),Abot(i,:))
    datamat3 = [datamat3, [[Xresc1(i,1:10:end)'; Xresc2(i,1:10:end)'], [ufull(i,1:10:end)'; flip(Abot(i,1:10:end))']]];
end
hold on 
plot(Xresc1(1,:),[u0steady; uf],'--')
csvwrite('uDiscreteTimestepsgamma20.csv',datamat3); 

% Compare with real steady state 

 