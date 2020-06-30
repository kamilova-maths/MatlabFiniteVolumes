function [ th, A, u, x, t ] = TimeDependentFDfull( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N )


% initialisation
th = zeros(N,K);
A  = zeros(N,K);
u  = zeros(N,K);

th(:,1) = th0;
A(:,1)  = A0;

% domain
dt = T/(K-1);
dx = L/N;
t = 0:dt:T;
x = ( dx:dx:L )';

% viscosity
mu = @(th) exp(gamma*th);

% heat source
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
 q = Q*(x>x1).*(x<x2);    % heat source
%q = zeros(N,1);

% laplacian for theta
Dx2 = spdiags( ones(N,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], N, N );
Dx2(N,N-1) = 2/dx^2;


% solve for u0
temp = 3*[D;A0].*mu([0;th0]);
tiph = ( temp(1:end-1) + temp(2:end) ) / 2;

Dx2u = spdiags( [ tiph(2:end), -(tiph(1:end-1)+tiph(2:end)), tiph(1:end-1) ] / dx^2, [-1,0,1], N-1, N-1 );
fu   = St*A0(1:N-1);
fu(1)  = fu(1)   - 1/D* tiph(1)  /dx^2;
fu(end)= fu(end) - 1*   tiph(end)/dx^2;

u(1:N-1,1) = Dx2u\fu;
u(N,1)     = 1;




for i=2:K
    %% Solve for A at next time step (explicit, hyperbolic)
    if sum( abs(u(:,i)) > dx/dt )
        err('CFL condition broken!')
    end
    
    tmpU = [1/D; u(:,i-1)];
    tmpA = [D; A(:,i-1);1];
    FL = tmpA(1:end-1).*tmpU;
    FR = tmpA(2:end  ).*tmpU;
    F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;
    F  = ( F(1:end-1) - F(2:end) ) / dx;
    A(:,i) = A(:,i-1) + dt * F;
    
    % possibly correct if A>1 somewhere (open can of worms)
        
        
    %% Solve for th at next time step (semi-implicit, parabolic)
    % Assemble matrices
    % - first order derivative
    temp = [ D; A(:,i); 1 ];                                
    g    = 1./A(:,i) .* ( temp(3:end) - temp(1:end-2) )/(2*dx) - u(:,i-1);  % notice u is treated explicitly
    Dx1 = spdiags( (ones(N,1).*g).*[ 1, -1 ]/(2*dx), [-1,1], N, N )';
    Dx1(N,N-1) = 0;
    % - reaction term
    Z = 2*Bi/Pe * spdiags( 1./sqrt(A(:,i)), 0, N, N );
    % - assemble matrix
    M = speye(N) - dt*( Dx1 + 1/Pe * Dx2 - Z );
    % - assemble rhs
     fth = (2*Bi/Pe ./sqrt(A(:,i-1))).*tha + q;
%    fth =  q;
    % - solve
    th(:,i) = M\( th(:,i-1) + dt*fth );
    
    
    
    %% Solve for u at next time step (laplacian)
    temp = 3*[D;A(:,i)].*mu([0;th(:,i)]);
    tiph = ( temp(1:end-1) + temp(2:end) ) / 2;

    Dx2u = spdiags( [ tiph(2:end), -(tiph(1:end-1)+tiph(2:end)), tiph(1:end-1) ] / dx^2, [-1,0,1], N-1, N-1 );
    fu   = St*A(1:N-1,i);
    fu(1)  = fu(1)   - 1/D* tiph(1)  /dx^2;
    fu(end)= fu(end) - 1*   tiph(end)/dx^2;

    u(1:N-1,i) = Dx2u\fu;
    u(N,i)     = 1;
    
end


th = [ zeros(1,K); ...
       th          ];

A  = [ D*ones(1,K); ...
       A           ];

u  = [ 1/D*ones(1,K); ...
       u           ];
   
x  = [0;x];
   
% to plot:
close all

subplot(1,3,1)
surf(t,x,th);
xlabel('$t$')
ylabel('$x$')
title('$\theta$')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

subplot(1,3,2)
surf(t,x,A);
xlabel('$t$')
ylabel('$x$')
title('$A$')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

subplot(1,3,3)
surf(t,x,u);
xlabel('$t$')
ylabel('$x$')
title('$u$')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

end
