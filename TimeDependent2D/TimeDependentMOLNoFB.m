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
N=4000; K=300;
% end of the domain
T = 1; L=1 ;

% We add the heaviside with H=1, and we remove it with H=0. 
H=0;
lam=1; 
% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;
% We try with tha=0
[P, A0, J0, th0, ~] = InitialConditionsSteady(K,gamma,Q,x1,x2,eps,St,tha,Bi,Pe,P0,R0,L,H,plt,lam);
 
%Initial conditions given by A0, th0
A0=A0';
th0=th0'; 
% Obtain initial condition for u
% viscosity
%mu = @(th) exp(-gamma*th); % Perhaps discretise this as well, but put a pin in that

dx = L/K;

uf= 1/A0(end);


y0(1:K) = A0;
%y0(1+K:2*K) = th0;
y0(1+K:2*K) = A0.*th0;


% Independent variable for ODE integration 
tout = linspace(0,T,N);

% ODE integration 
reltol = 1.0e-04; abstol = 1.0e-04;
options = odeset('RelTol',reltol,'AbsTol',abstol);
%[t,y] = ode15s(@coupledPdeNoFB,tout,y0);
[t,y] = ode15s(@coupledPdeNoFBwithFV,tout,y0);

A = y(1:N,1:K); % N rows, K columns, each row is a timestep, each column is space point
w = y(1:N,K+1:2*K);

th = w./A; 
%% Solve for u at next time step (laplacian)
%     temp = 3*mu([0;0;th(:,i)]);     % add ghost node to th
%     tiph = ( temp(1:end-1) + temp(2:end) ) / 2;
u = zeros(size(A));
u0 = usolutionNoFB(A0,th0);  % Do I even need this guy? 
u(1,:) = u0; 
for i=2:N
    u(i,:) = usolutionNoFB(A(i,:)',th(i,:)');   
    
end

th = [ zeros(N,1), ...
       th          ];

A  = [ D*ones(N,1), ...
       A           ];
   
u  = [ u , ...
       uf.*ones(N,1) ];
   
dx = 1/K;

x = (0:dx:L)';

% 
figure;
surf(t,x,A','LineStyle','none')
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
title('$A$ with MOL','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% figure;
% surf(t,x,thfull')
% xlabel('$t$')
% ylabel('$x$')
% title('$\theta$ with MOL')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% figure;
% surf(t,x,ufull')
% xlabel('$t$')
% ylabel('$x$')
% title('$u$ with MOL')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)

% Compare with previous code
plt = 0; 
th0FD = th(end,2:end)';
A0FD  = A(end,2:end)'; 


[ th1, A1, u1, x1, t1 ] = TimeDependentFDfull_v3( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N, plt);

% figure;
% surf(t,x,A1)
% xlabel('$t$')
% ylabel('$x$')
% title('$A$ with explicit discretisation')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% figure;
% surf(t,x,th1)
% xlabel('$t$')
% ylabel('$x$')
% title('$\theta$ with explicit discretisation')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% figure;
% surf(t,x,u1)
% xlabel('$t$')
% ylabel('$x$')
% title('$u$ with explicit discretisation')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)

figure; plot(x,A(end,:))
hold on 
plot(x,A1(:,end),'--','LineWidth',2)
legend({'MOL','Explicit FD'},'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

figure; plot(x,th(end,:))
hold on 
plot(x,th1(:,end),'--','LineWidth',2)
legend({'MOL','Explicit FD'},'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

figure; plot(x,u(end,:))
hold on 
plot(x,u1(:,end),'--','LineWidth',2)
legend({'MOL','Explicit FD'},'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)