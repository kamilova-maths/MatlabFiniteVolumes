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
%P0 = (10000*L)/((R1^2)*uc*mu0);
%P0=0.9427;

DeltaT = (Qc*L)/(rho*c*uc);
Bi= ((L^2)*h)/(k*R1); 
tha = 0.005; 
D = (R0^2)/(R1^2); 

dudx0=-0.7; 
P0 =- 3*(1)*D*dudx0; 


gamma = 0; 

%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
eps = 1e-2;

% Calculating the initial conditions as a solution of the steady state
% problem 
K=2500; N=100;
% end of the domain
T = 5; L=1; 


[P, A0, J0, th0] = InitialConditionsSteady(N,gamma,Q,x1,x2,eps,St,tha,Bi,Pe,P0,R0);

% we change the dimensions so that it is compatible with out FD code. 
A0=A0';
th0=th0'; 
xplot=linspace(0,1,N); 
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

% tests for boundary conditions
dx=L/N;

dAdx=derivative(A0,dx); 

lhs = P'./(3.*exp(gamma.*th0));

rhs = (1./A0).*dAdx';
abs(lhs(1)-rhs(1))
