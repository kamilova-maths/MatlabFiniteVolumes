function [ lambda, A, u, x, t ] = TimeDependentFDwithLambdatConstantmu( A0, D, P0, St, T, L, K, N, plt)

% initialisation
%th = zeros(N,K);
A  = zeros(K,N);
u  = zeros(K,N);
lambda = zeros(N,1); 

%th(:,1) = th0;
A(:,1)  = A0;

% domain
dt = T/(N-1);
dx = L/K;
t = 0:dt:T;
x = ( dx:dx:L )'; % Since we have Dirichlet boundary conditions, we don't need x=0
I = find(A0>0.999,1,'first');
lambda(1)= x(I);  % this is an initial guess. If it doesn't work, it won't be because of this .

% viscosity - constant, in this case 
%mu = @(th) exp(-gamma*0); % always going to be 1 in this version 
mu=1; 
% heat source

 
%q = zeros(K,1);

% solve for u0

A0  = [ D; A0  ];           % add ghost node to A


uf= 1/A0(end);


Dx2u = spdiags( [ A0(2:end).*3*mu, -(A0(1:end-1).*3*mu+A0(2:end).*3*mu), A0(1:end-1).*3*mu ] / dx^2, [-1,0,1], K, K );
Dx2u(1,2) = Dx2u(1,2) + A0(1).*3*mu / dx^2;        % include effect from Neumann BC
fu   = - St*( A0(1:end-1) + A0(2:end) )/ 2;
%fu(1)  = fu(1)   - 2*P0/(3*D*dx); % include derivative (again, Neumann BC)
fu(1) = fu(1) -2*3*mu*P0(1)*lambda(1)/(3*dx); 

%fu(1) =  fu(1) - 2*(1/A0(1)*(A0(3)-A0(2))/(dx)); 
%fu(end)= fu(end) - 1   * A0(end)* tiph(end)/dx^2;

fu(end)= fu(end) - uf   * A0(end)* 3*mu/dx^2;

u(:,1) = Dx2u\fu;

dudx = derivative(u(:,1),dx)'; 
dAdx = derivative(A(:,1),dx)';
lambda(2) = dt.*(u(I,1) + (1.*(dudx(I)./dAdx(I)))) + lambda(1);
ldot =  1+ (u(I,1) - u(I-1,1))/(A(I,1)-A(I-1,1));
for i=2:N
    % Change these to value forms
    
    %% Solve for A at next time step (explicit, hyperbolic)
    if sum( abs(u(:,i-1)) > dx/dt )
        error('CFL condition broken!')
    end

    tmpU = [ u(:,i-1); uf ] -ldot.*[x; 1] ; % size K+1 x 1 
    tmpA = [D; A(:,i-1);1];
    FL = tmpA(1:end-1).*tmpU;
    FR = tmpA(2:end  ).*tmpU;
    % impose boundary conditions correctly ... 
    F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;
    F  = ( F(1:end-1) - F(2:end) ) / dx; 
    % Source term
   
    
    S  = -A(:,i-1).*ldot./lambda(i); 
    A(:,i) = A(:,i-1) + dt*((1/lambda(i)).*F+S);
    % possibly correct if A>1 somewhere (open can of worms)
              
    ldot = 1+ (u(K,i) - u(K-1,i))/(A(K,i)-A(K-1,i));
    %% Solve for u at next time step (laplacian)
%     temp = 3*mu([0;0;th(:,i)]);     % add ghost node to th
%     tiph = ( temp(1:end-1) + temp(2:end) ) / 2;
    
    % Technically can do this in one line, but I can't figure out how -
    % this gives same results as above so we don't bother changing it
   
    A0  = [ D; A(:,i)  ];           % add ghost node to A
    %A0 =  [A(:,i); 1]; 
    Dx2u = spdiags( [ A0(2:end).*3.*mu, -(A0(1:end-1).*3.*mu+A0(2:end).*3*mu), A0(1:end-1).*3.*mu ] / dx^2, [-1,0,1], K, K );
    Dx2u(1,2) = Dx2u(1,2) + A0(1).*3.*mu / dx^2;              % include effect from Neumann BC
    fu   = - (lambda(i).^2)*St*( A0(1:end-1) + A0(2:end) )/ 2;

    if length(P0)==1
    fu(1) = fu(1) -2*3*mu*P0*lambda(i)/(3*dx);  % include derivative (again, Neumann BC)
    else
    fu(1) = fu(1) -2*3*mu*P0(i)*lambda(i)/(3*dx);  % include derivative (again, Neumann BC)
    end
    
    fu(end)= fu(end) - uf   * A0(end)* (3*mu)/dx^2;
    u(:,i) = Dx2u\fu;
    
    dudx = derivative(u(:,i),dx)'; 
    dAdx = derivative(A0,dx)';
   
    lambda(i+1) = dt.*(u(end,i) + (A0(end).*(dudx(end)./dAdx(end)))) + lambda(i);
end


A  = [ D*ones(1,N); ...
       A           ];

% u  = [ 1/D*ones(1,K); ...
%        u           ];
u  = [ u ; ...
       uf.*ones(1,N)   ];

   
% u = [u; ...
%      ufinal.*ones(1,K)];
   
x  = [0;x];% We add zero at the end
   
% to plot:
close all


if plt == 1

subplot(1,2,1)
surf(t,x,A);
xlabel('$t$')
ylabel('$x$')
title('$A$')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

subplot(1,2,2)
surf(t,x,u);
xlabel('$t$')
ylabel('$x$')
title('$u$')
set(gca,'TickLabelInterpreter','latex','fontsize',13)
end

end

