% Initialise all the values I need to run time dependent finite difference
% code 

% Here we set up all the parameters, and we also calculate the initial
% conditions as a solution of the steady state problem 


close all 
clear all 

clc

set(0,'DefaultAxesFontSize',12,'DefaultTextInterpreter','latex');

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
P0val = (10000*Ld)/((R1^2)*uc*mu0);

DeltaT = (Qc*Ld)/(rho*c*uc);
Bi= ((Ld^2)*h)/(k*R1); 
tha = 0.005; 
D = (R0^2)/(R1^2); 

gamma = 30; 



% Calculating the initial conditions as a solution of the steady state
% problem 
K=4000; N=100;
% end of the domain
T = 1; L=1 ;

%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Qvalue = 1;
x=linspace(0,L,N);
Q = Qvalue*(x>x1).*(x<x2);    % heat source  
eps = 1e-2;

% We add the heaviside with H=1, and we remove it with H=0. 
H=0;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;

%[P, A0, J0, th0,uf] = InitialConditionsSteady(N,gamma,Qvalue,x1,x2,eps,St,tha,Bi,Pe,P0val,R0,L,H,plt);

% we change the dimensions so that it is compatible with out FD code. 

% Calculating the initial conditions as a solution of the steady state
% problem 
K=8000; N=200;
% end of the domain
T = 1; L=1 ;

% We add the heaviside with H=1, and we remove it with H=0. 
H=0;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 1;
% We try with tha=0
[P, A0, J0, th0,uf] = InitialConditionsSteady(N,gamma,Qvalue,x1,x2,eps,St,tha,Bi,Pe,P0val,R0,L,H,plt);

% we change the dimensions so that it is compatible with out FD code. 
x=linspace(0,L,N);
t = linspace(0,T,K);

A0=A0';
th0=th0'; 
u0 = (1./A0);
P0 = P0val*ones(1,K); 
[ th1, A1, u1, x1, t1 ] = TimeDependentFDfull_v3( th0, A0, A0(1), gamma, P0, Pe, St, Bi, tha, T, L, K, N, uf, plt );
A0 = A1(2:end,end);
th0= th1(2:end,end);
omega = 2*pi; 
%P0 = P0val*exp(omega.*t); 
P0 = P0val+ P0val*sin(omega.*t); 
[ th2, A2, u2, x2, t2 ] = TimeDependentFDfull_v3( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N, uf, plt );

A = A2(2:end,:);
th = th2(2:end,:);
u = u2(2:end,:);

return
%Export plotted data

valuesmatrix=[x',A(:,1), A(:,1600:1600:end)];

csvwrite('AwithP0tomega2pi.csv',valuesmatrix)

valuesmatrix=[x',u(:,1), u(:,1600:1600:end)];

csvwrite('uwithP0tomega2pi.csv',valuesmatrix)

valuesmatrix=[x',th(:,1), th(:,1600:1600:end)];

csvwrite('thwithP0tomega2pi.csv',valuesmatrix)

%dlmwrite('thwithP0tomega4.csv', valuesmatrix,'delimiter',',','precision',4)

valuesP0 = [t(1:10:end)', P0(1:10:end)']; 
csvwrite('P0withomega2pi.csv',valuesP0)