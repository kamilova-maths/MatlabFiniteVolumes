% This code is very similar to InitialiseforFD, but here we allow to vary
% one of the parameters. The steps will be the same. 


% General parameters definition

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
Qvalue = 1;
%Q = Qvalue*(x>x1).*(x<x2);    % heat source  
eps = 1e-2;

% Calculating the initial conditions as a solution of the steady state
% problem 
K1=1000; N=100;
% end of the domain
T = 1; L=1 ;

% Control parameters that set plots and use of Heaviside

% We add the heaviside with H=1, and we remove it with H=0. 
H=0;

% Plots for steady state - 1 , no plots for steady state - 0
pltSt = 0;
pltTD = 0; 
% Parameter values we will loop over

numeigens = 20; 
% For Peclet
%parVal = [10, 20, 30, 40, 50];
% For Bi
%parVal = [10, 50, 100, 150, 200];
% For St
parVal = [20,25,30,32];
eigVals = zeros(numeigens,length(parVal)); 
for i=1:length(parVal)

%Pe = parVal(i); 
%Bi = parVal(i);
St = parVal(i);
%obtaining the steady state, with K1 (smaller than K2) 
[P, A0, J0, th0,uf] = InitialConditionsSteady(N,gamma,Qvalue,x1,x2,eps,St,tha,Bi,Pe,P0,R0,L,H,pltSt);

A0=A0';
th0=th0'; 
u0 = (1./A0);

[ th1, A1, u1, ~, t1 ] = TimeDependentFDfull_v3( th0, A0, A0(1), gamma, P0, Pe, St, Bi, tha, T, L, K1, N, uf, pltTD );

%This gives us the steady state we now use to compute the real thing

A0 = A1(2:end,end);
th0= th1(2:end,end); 
K2 = 3*K1; 

[ th2, A2, u2, ~, t2 ] = TimeDependentFDfull_v3( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K2, N, uf, pltTD );

th0 = th2(:,1);
A0 = A2(:,1);
u0 = u2(:,1);

pltEig=1;

eigVals(:,i) = LinearStabilityFun(th0,A0,u0,gamma,Pe,Bi,St,L,x1,x2,pltEig,tha,numeigens) ;
end

plot(eigVals,'o')
%legend({'Pe$= 10$','Pe$= 20$','Pe$= 30$','Pe$= 40$','Pe$= 50$'},'Interpreter','latex')
%legend({'Bi$= 10$','Bi$= 50$','Bi$= 100$','Bi$= 150$','Bi$= 200$'},'Interpreter','latex')
%legend({'St$= 1$','St$= 5$','St$= 10$','St$= 15$','St$= 20$','St$= 25$'},'Interpreter','latex')
legend({'St$= 20$','St$= 25$','St$= 30$','St$= 32$'},'Interpreter','latex')



set(gca,'TickLabelInterpreter','latex','fontsize',15)
xlabel('Re($\sigma$)')
ylabel('Im($\sigma$)')