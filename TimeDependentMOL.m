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
Ta = 0.005; 
D = (R0^2)/(R1^2); 

alpha = 50; 

%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
p.eps = 1e-2;
p1 = p;


% number of iterations
n=20; 

% Initial conditions - one for each variable 

[Pinit, Ainit, Jinit,Thinit]=InitialConditionsSteady(n);





% Independent variable for ODE integration
t0=0.0;
tf = 3; 
tout = linspace(t0,tf,n);

nout=n;


% ODE integration 
reltol = 1.0e-08; abstol=1e-08; 
options = odeset('RelTol',reltol,'AbsTol',abstol); 

% My dimensions are wrong, I need to think about this a bit more. Maybe I
% will do them separately? 

% for A

N = 2000;
K= 2100; 
x=linspace(0,1,N);
t=linspace(0,1,K);
dx = 1/(N-1);
dt = 1/(K-1); 

% WE WILL HAVE A MATRIX WITH ROWS CORRESPONDING TO EACH T, AND COLUMNS
% CORRESPONDING TO EACH X. THIS WILL NOT CHANGE. 
A = zeros(K,N);
P = zeros(K,N);
J = zeros(K,N);
u = zeros(K,N); 
Th = zeros(K,N); 
mu = ones(1,N);

Qval=1;

Q = Qval*(x>x1).*(x<x2);    % heat source

thetaa=Ta*ones(1,N);
 % For the first timestep, we use the steady state solutions.
 % Essentially, we are using them as initial conditions
[Pold, Aold, Jold,Thold]=InitialConditionsSteady(N);
uold = 1./Aold;
% We save the initial conditions in k=1, for all n, for each solution
% matrix. 

A(1,:)=Aold;
P(1,:)=Pold;
u(1,:)=uold;
Th(1,:)=Thold; 
k=1; n=1;  

% We have decided to split the problem into two parts. The first part is
% when A<1. Once it reaches this value (or surpasses it), we switch to the
% second problem where we impose A to be 1. 
while(A(k,n)<1)
    
for k=2:K % We loop through each timestep


 for n=1:N % We loop though each x value - for a fixed timestep
     % We setup the boundary conditions 
        if n==1
            A(k,n)= D; 
            u(k,n) = 1/D;
            P(k,n)=P0; 
            Th(k,n)= 0;
        % We set the N+1 imaginary gridpoint to zero    
        elseif n==N
            A(k,n)=(dt/3*mu(n))*P(k-1,n)+((dt/dx)*u(k-1,n)+1)*A(k-1,n); 
            Th(k,n) = (dt/Pe)*(-1/dx + (2*dx+Th(k-1,n-1))/dx^2)+(2*Bi*dt/Pe)*(1/sqrt(A(k-1,n)))*(dx+thetaa(n))+(Q(n)-u(k-1,n))*dt-dx; 
            P(k,n) = St*A(k,n-1)*dx +P(k,n-1); 
            u(k,n)=-dx*P(k,n-1)/(3*mu(n-1)*A(k,n-1))+u(k,n-1); 
        else
        % This populates the matrix for A and theta. This gives us A(x,t) at each
        % timestep, so at the end of the for loop we have all the
        % timesteps.
        
        % By the time you get to this part of the code, k is larger or
        % equal to 2, and n is larger or equal to 2 but less than N 
       
         A(k,n) =  (dt/3*mu(n))*P(k-1,n) +((dt/dx)*u(k-1,n)+1)*A(k-1,n)-(dt/dx)*u(k-1,n)*A(k-1,n+1);

         dadx=(A(k-1,n+1)-A(k-1,n))/dx;
         dthdx=(Th(k-1,n+1)-Th(k-1,n))/dx;
         ddthddx = (Th(k-1,n+1)-2*Th(k-1,n)+Th(k-1,n-1))/dx^2; 
         Th(k,n) = dt/Pe*((1/A(k-1,n))*dadx*dthdx+ddthddx) - (2*Bi*dt/Pe)*(1/sqrt(A(k-1,n)))*(Th(k-1,n)-thetaa(n))+...
             (Q(n)-u(k-1,n)*dthdx)*dt+Th(k-1,n);
         % I think these following two are correct, but that would mean
         % that all the 'n' above are actually n-1. I am looping over the
         % 'n', so they are increasing. In this section, they are between 2
         % and N-1. So there might be a different rule for P and u, than
         % there is for A and Theta - I have to think about it some more
         
         % P of the next n is going to depend on P of the previous n and A
         % of the previous n, all at the same k. 
         
         %He = max(1-exp(-(1-A(k,n))/eps),0);% regularised Heaviside function
        % P(k,n)=St*A(k,n-1)*dx*He+P(k,n-1);
         P(k,n)=St*A(k,n-1)*dx+P(k,n-1);
         
         u(k,n)=-dx*P(k,n-1)/(3*mu(n-1)*A(k,n-1))+u(k,n-1); 
       
         % if A reaches one, we will change the problem to solve for the
         % rest of the domain. We are assuming here that once it reaches
         % one it doesn't go back to being less than one. 
        
        end
end
     
     
end
         
         % PLEASE MAKE A STOP AT 1. IT CURRENTLY DOES NOT STOP AT 1 AND
         % EVERYTHING BLOWS UP. 
       
end

% Remember that if we broke the while loop, we have reached A larger than
% one. Therefore, we have to replace A(n,k) with 1. 

% First we have to finish the current loop, then we continue with the
% normal loops. 
% But this means that we also have to chop the domain. Maybe this isn't
% good enough. 
return


          if A(k,n)>1
            prb=2; 
         end
         
         if prb==2
        A(k,n)=1;
        P(k,n)=St*dx+P(k,n-1);
        u(k,n)=-dx*P(k,n-1)/(3*mu(n-1))+u(k,n-1); 
        ddthddx = (Th(k-1,n+1)-2*Th(k-1,n)+Th(k-1,n-1))/dx^2; 
        Th(k,n) = (dt/Pe)*(ddthddx) - (2*Bi*dt/Pe)*(Th(k-1,n)-thetaa(n))+...
             (Q(n)-u(k-1,n)*dthdx)*dt+Th(k-1,n); 
     end

[T,X]=meshgrid(t,x);
contour(T,X,A'); colorbar;
set(gca,'YDir','reverse')
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
% [tA,A] = ode15s(@(t,A) SolvingODEforA(t,A,n,mu,u,P),tout,Ainit,options);
% 
% % Now we compute the next timestep for P and u, using the solution for A
% 
% 
% % for J
% opt=3;
% [tJ,J] = ode15s(@(t,y) SteadyStateMOL(t,y,n,opt),tout,Jinit,options);
% % for Theta
% opt=4;
% [tTh,Th] = ode15s(@(t,y) SteadyStateMOL(t,y,n,opt),tout,Thinit,options);
% 
% % The results above are of the format: 
% % ROWS are values at particular t, 
% % COLUMNS are values at particular x
% 
% % In order to plot, you must plot the transverse of each matrix above. 
% 
% tplot=linspace(0,1,n);
% 
% subplot(2,2,1)
% plot(tplot,P')
% xlabel('$t$','Interpreter','latex', 'Fontsize', 15)
% ylabel('$P(t)$','Interpreter','latex','Fontsize',15)
% set(gca,'TickLabelInterpreter','latex','fontsize',15)
% 
% subplot(2,2,2)
% plot(tplot,A')
% xlabel('$t$','Interpreter','latex', 'Fontsize', 15)
% ylabel('$A(t)$','Interpreter','latex','Fontsize',15)
% set(gca,'TickLabelInterpreter','latex','fontsize',15)
% 
% subplot(2,2,3)
% plot(tplot,J')
% xlabel('$t$','Interpreter','latex', 'Fontsize', 15)
% ylabel('$J(t)$','Interpreter','latex','Fontsize',15)
% set(gca,'TickLabelInterpreter','latex','fontsize',15)
% 
% 
% subplot(2,2,4)
% plot(tplot,Th')
% xlabel('$t$','Interpreter','latex', 'Fontsize', 15)
% ylabel('$\theta(t)$','Interpreter','latex','Fontsize',15)
% set(gca,'TickLabelInterpreter','latex','fontsize',15)
% 
% 
% 
% 
% 
