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
K=4000; N=100;
% end of the domain
T = 1; L=1 ;

% We add the heaviside with H=1, and we remove it with H=0. 
H=0;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;
% We try with tha=0
[P, A0, J0, th0,uf] = InitialConditionsSteady(N,gamma,Q,x1,x2,eps,St,tha,Bi,Pe,P0,R0,L,H,plt);

% we change the dimensions so that it is compatible with out FD code. 
x=linspace(0,L,N);

A0=A0';
th0=th0'; 
u0 = (1./A0);

plt = 0; 
[ th1, A1, u1, x1, t1 ] = TimeDependentFDfull_v3( th0, A0, A0(1), gamma, P0, Pe, St, Bi, tha, T, L, K, N, uf, plt);

Ast = A1(1:end,end);
tst = th1(1:end,end);
ust = u1(1:end,end);

%A0  = A1(2:end,end)+[0; sin(x(1:end-2)');0];
A0  = A1(2:end,end);

%th0 = th1(2:end,end)+[0; 0.1*exp(-x(1:end-2)');0];
th0 = th1(2:end,end);
u0  = u1(2:end,end);
plt = 1;
[ th2, A2, u2, x2, t2 ] = TimeDependentFDfull_v3( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N, uf, plt);

return

% for plotting 
figure; 
dx=1/N;
xplot=0:dx:L;
xplot=xplot';


for i=1:100:K
   if i==1
   plot(xplot,A2(:,i),'--','LineWidth',2)    
   else
   plot(xplot,A2(:,i),'HandleVisibility','off')
   end
   hold on     
end
plot(xplot,Ast,'--','LineWidth',2)
xlabel('$x$')
ylabel('$A$')
legend({'Initial condition','Steady state'},'Interpreter','latex')
title('A each 50 timesteps')
% 
 set(gca,'TickLabelInterpreter','latex','fontsize',15)
 
 figure; 
 for i=1:100:K
   if i==1
   plot(xplot,u2(:,i),'--','LineWidth',2)    
   else
   plot(xplot,u2(:,i),'HandleVisibility','off')
   end
   hold on     
end
plot(xplot,ust,'--','LineWidth',2)
xlabel('$x$')
ylabel('$u$')
legend({'Initial condition','Steady state'},'Interpreter','latex')
title('u each 50 timesteps')
% 
 set(gca,'TickLabelInterpreter','latex','fontsize',15)

 

 figure; 
 for i=1:100:K
   if i==1
   plot(xplot,th2(:,i),'--','LineWidth',2)    
   else
   plot(xplot,th2(:,i),'HandleVisibility','off')
   end
   hold on     
end
plot(xplot,tst,'--','LineWidth',2)
xlabel('$x$')
ylabel('$\theta$')
legend({'Initial condition','Steady state'},'Interpreter','latex')
title('$\theta$ each 50 timesteps')
% 
 set(gca,'TickLabelInterpreter','latex','fontsize',15)


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

% surface plot - so that you don't have to rerun the whole thing
dt = T/(K-1);
dx = L/N;
t = 0:dt:T;
x = ( 0:dx:L )'; % Since we have Dirichlet boundary conditions, we don't need x=0 
figure(1)
subplot(1,3,1);
h1=surf(t,x,th2);
set(h1,'LineStyle','none')
xlabel('$t$')
ylabel('$x$')
title('$\theta$')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

subplot(1,3,2);
h2=surf(t,x,A2);
set(h2,'LineStyle','none')
xlabel('$t$')
ylabel('$x$')
title('$A$')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

subplot(1,3,3);
h3=surf(t,x,u2);
set(h3,'LineStyle','none')
xlabel('$t$')
ylabel('$x$')
title('$u$')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

print(gcf,'-dpdf','test1.pdf','-bestfit')