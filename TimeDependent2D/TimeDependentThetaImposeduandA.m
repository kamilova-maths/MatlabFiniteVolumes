% Code to compute time dependent equation for theta, where we impose u and
% A, in order to find the right discretisation procedure for the equations
% that we have. 

% This is a MOL code, where theta is time dependent and spatial dependent. 
% We have a transport term, that we anticipate will give us problems, a
% diffusion term that should be fine, provided the Peclet number isn't
% huge compared to the other terms, a reaction term which is just
% multiplied by an identity and an easy source term which will make up the
% right hand side. 

% First we compute the steady state solution, which gives us a value for
% the known parameters in the problem (and initial value) 

% We define the dimensional parameters
close all 
clear all 

clc

set(0,'DefaultAxesFontSize',12,'DefaultTextInterpreter','latex');

rho= 1.8*10^3; %Bergstrom ; 
g = 10; 
c= 900; % Fitt and Howell
L=7; % Temperature Profiles in Soderberg Electrodes
u = 10^-5; %Bergstrom approximation
R0=0.5;
R1=1; 
k = 3; 
Qc=15000; %Taken very vaguely from Temperature profiles in Soderberg electrodes. 
mu0=10^10; % given by Bjornar at a reference temperature
T_a=343; 
h = 7; 

% We define the non-dimensional parameters

Pe = (rho*c*u*L)/(k);
%Pe = 80;
epsilon=R1/L;
St=(rho*g*L^2)/(u*mu0);
P0 = (10000*L)/((R1^2)*u*mu0);

Bi= ((L^2)*h)/(k*R1); 

DeltaT = (Qc*L)/(rho*c*u);

% We define parameters for the problem, source term and convective heating

x1 = 5/7;
x2 = 6.5/7;
Q = 1;
eps = 1e-2;

% temperature dependence for viscosity 
alpha = 50; 

% constant heating at the wall

Ta = 0.001;

% we define the values for N  - discretisation in SPACE, and K-
% discretisation in TIME

N = 1000;
K= 1500; 
x=linspace(0,1,N);
t=linspace(0,1,K);

% I am not sure if dt and dx are quite right 
dx = 1/(N+1);
dt = 1/(K+1); 

% We obtain the initial conditions, which will be given as vectors of size
% 1xK 
[PinitK, AinitK, JinitK, ThinitK]=InitialConditionsSteady(K,alpha,Q,x1,x2,eps,St,Ta,Bi,Pe,P0,R0);

% We obtain the initial conditions, which will be given as vectors of size
% 1xN
[PinitN, AinitN, JinitN, ThinitN]=InitialConditionsSteady(N,alpha,Q,x1,x2,eps,St,Ta,Bi,Pe,P0,R0);


% We calculate the values of g and f, according to notes. In particular, 
% they could be something way easier, like a straight line, but I think
% that this will help me implement the actual values later on when coupling
% everything together.

% Both g and f must be vectors of size 1xN

g = (1./AinitN).*((1/Pe)*derivative(AinitN,dx)-1); 
% g is negative. Just pointing it out. 

Qfun = Q*(x>x1).*(x<x2);    % heat source


% Now we write the discretisation matrices, based on soltheta.m in
% HeatEquationMOL folder. 
% second order x diffusion operator, homogenous Dirichlet at x = 0,
% homogenous Neumann at x=1
Dx2 = sparse(1:N,1:N,-2/dx^2*ones(N,1),N,N) ...
    + sparse(1:N-1,2:N,1/dx^2*ones(N-1,1),N,N) ...
    + sparse(2:N,1:N-1,1/dx^2*[ones(N-2,1); 2],N,N); 
% Dx2 = spdiags( ones(N,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], N, N );
% Dx2(N,N-1) = 2/dx^2;

Dx2(1,:)=0;

M2 = kron(speye(K),Dx2);

% Transport term 

Dx1 = sparse(1:N,1:N,1/dx*ones(N,1),N,N) ...
    + sparse(2:N,1:N-1,-1/dx*ones(N-1,1),N,N); % forward difference x derivative, homogeneous Dirichlet at x = 0
Dx1(end,:)=0; 


gDx1 = Dx1.*g'; 
gM3 = kron(speye(K),gDx1);


% time derivative
Dt1 = sparse(1:K,1:K,1/dt*ones(K,1),K,K) ...
    + sparse(2:K,1:K-1,-1/dt*ones(K-1,1),K,K);


M1 = kron(Dt1,speye(N)); 
z= (2*Bi/Pe).*(1./AinitN.^(1/2)); 
Imat=kron(speye(K),speye(N).*z');

% We want to satisfy the boundary condition, so theta has to be zero when
% x=0. 

% Now we can build the matrix that discretises the whole lhs of the
% equation
% indices of nodes - helpful way to write them out 

is = 1:N*K; % list of all the nodes
ist0 = 1:N; % nodes corresponding to t=0; 
isx0 = 1:N:N*K; % nodes corresponding to x=0;
isy1 = N:N:N*K; % nodes corresponing to x=1. 

% Discretised version of PDE M*theta=b

M=M1-gM3-(1/Pe)*M2+Imat;

% We want to satisfy the boundary condition, so theta has to be zero when
% x=0. 
f = (2*Bi/Pe).*(1./(AinitN.^(1/2)))*Ta+Qfun; 
%f(1)=0; 
b=repmat(f,[1,K])'; 

% initial and boundary conditions 

% values for t = 0
thetat0 = AinitN; 

% values for x = 0
thetax0 = zeros(1,K); 


% solve equation 
is_prescribed = [ist0 isx0]; %nodes at which theta is prescribed
is_solve1 = setdiff(is,is_prescribed); % nodes at which to solve for theta
is_solve2 = setdiff(is,ist0); % rows of matrix to use
theta = NaN*ones(N*K,1); % initialise 
theta(ist0)=thetat0; % initial conditions 
theta(isx0)=thetax0; % boundary condition 
theta(is_solve1)=M(is_solve2,is_solve1)\(b(is_solve2) - M(is_solve2,is_prescribed)*theta(is_prescribed)); 


theta_calc=reshape(theta,N,K); 
theta_calc=theta_calc';
figure(1)
mesh(x(:),t(:),theta_calc)

ylabel('$t$'); xlabel('$x$'); 

set(0,'DefaultAxesFontSize',12,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14)

