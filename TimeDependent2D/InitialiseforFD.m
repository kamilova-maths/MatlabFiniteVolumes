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

gamma = 30; 

%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
eps = 1e-2;

% Calculating the initial conditions as a solution of the steady state
% problem 
K=4000; N=300;
% end of the domain
T = 1; L=1 ;

% We add the heaviside with H=1, and we remove it with H=0. 
H=0;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;
% We try with tha=0
[P, A0, J0, th0,uf] = InitialConditionsSteady(N,gamma,Q,x1,x2,eps,St,tha,Bi,Pe,P0,R0,L,H,plt);

% we change the dimensions so that it is compatible with out FD code. 
A0=A0';
th0=th0'; 
xplot=linspace(0,1,N); 


[ th1, A1, u1, x1, t1 ] = TimeDependentFDfull_v3( th0, A0, A0(1), gamma, P0, Pe, St, Bi, tha, T, L, K, N, uf );

A0 = A1(2:end,end);
th0= th1(2:end,end); 


[ th2, A2, u2, x2, t2 ] = TimeDependentFDfull_v3( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N, uf );



% for plotting 
% dx=1/N;
% xplot=0:dx:L;
% xplot=xplot';
% 
% 
% for i=1:50:K
%     
%    plot(xplot,A(:,i))
%    hold on 
%    
%     
% end
% 
% xlabel('$x$')
% ylabel('$A$')
% title('A each 50 timesteps')
% 
% set(gca,'TickLabelInterpreter','latex','fontsize',15)

% % tests for boundary conditions
% dx=L/N;
% 
% dAdx=derivative(A0,dx); 
% 
% lhs = P'./(3.*exp(gamma.*th0));
% 
% % rhs = (1./A0).*dAdx';
% % abs(lhs(1)-rhs(1));
% close all 
% figure(1)
% for i=1:5:K
%    plot(x1,th1(:,i))
%    hold on
%    pause
% end 
