function [err,theta_ex,theta_calc] = MSTDThetaImposeduandA(n,Pe)

% Manufactured solutions - ONLY WORKS FOR N=K, to not struggle with x vs 
%t discretisations: 
% Here we are trying to replicate the simplest possible case, which is when
% theta is just a function of x, and g,z are just one. This way, f is also
% just a function of x (and we don't have to play with matrices), and our
% initial condition is just x. We take N=K so that the initial condition
% and boundary conditions are compatible. The problem we are trying to
% recover here is

% theta= x^2-2x, with f=x^2-4x, theta=0 at x=0 and theta=x^2-4x at t=0, as well as
% dtheta/dx=0 at x = 1. 

% We are basing the discretisation procedure on soltheta.m 
close all 
% Easy Ainit


N=n;
K=N;


x=linspace(0,1,N);
t=linspace(0,1,K);

% This first run is for the simplest case of g=z=1, so we just take them
% out of the equation altogether


% I am not sure if dt and dx are quite right 
dx = 1/(N+1);
dt = 1/(K+1); 


Dx2 = sparse(1:N,1:N,-2/dx^2*ones(N,1),N,N) ...
    + sparse(1:N-1,2:N,1/dx^2*ones(N-1,1),N,N) ...
    + sparse(2:N,1:N-1,1/dx^2*[ones(N-2,1); 2],N,N); 

% boundary condition at x=0. The condition at x=1 is still homogeneous
% Neumann
Dx2(1,:)=0;

M2 = kron(speye(K),Dx2);

Dx1 = sparse(1:N,1:N,1/dx*ones(N,1),N,N) ...
    + sparse(2:N,1:N-1,-1/dx*ones(N-1,1),N,N); % forward difference x derivative, homogeneous Dirichlet at x = 0
% Homogeneous Neumann condition at x=1
Dx1(end,:)=0; 


M3 = kron(speye(K),Dx1); 


% time derivative
Dt1 = sparse(1:K,1:K,1/dt*ones(K,1),K,K) ...
    + sparse(2:K,1:K-1,-1/dt*ones(K-1,1),K,K);
% time derivative with Homogeneou sDirichlet at t=0 

M1 = kron(Dt1,speye(N)); 


Imat=kron(speye(K),speye(N));

% I somehow feel that this is still right, otherwise I am imposing theta =
% 1 in places where I know I don't need to solve. 
% An alternative would be to change this in M, provided I am SURE which
% values are which in M. 

% indices of nodes - helpful way to write them out 

is = 1:N*K; % list of all the nodes
ist0 = 1:N; % nodes corresponding to t=0; 
isx0 = 1:N:N*K; % nodes corresponding to x=0;
isy1 = N:N:N*K; % nodes corresponing to x=1. 

% Discretised version of PDE M*theta=b

M=M1-M3-(1/Pe)*M2+Imat;

% We make matrix of f, and for the case f only depends on one variable, set
% the other to zero. It is easier to have it all imposed this way. 
xmat = repmat(x,K,1)'; 

tmat = repmat(t,N,1); 

fmat = 1-4*xmat+xmat.^2+tmat;

b=fmat(:);

% initial and boundary conditions 

% values for t = 0
thetat0 = x.^2-2*x;

% values for x = 0
thetax0 = t; 

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





figure(2)
theta_ex=tmat+xmat.^2-2*xmat;
mesh(x(:),t(:),theta_ex')


err=0; 




end