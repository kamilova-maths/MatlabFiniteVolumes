function [ th, A, u, x, t ] = TimeDependentFDfull_v4( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N )


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
mu = @(th) exp(-gamma*th);

% heat source
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
 q = Q*(x>x1).*(x<x2);    % heat source  % CHANGE ME BACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%q = zeros(N,1);


% laplacian for theta
Dx2 = spdiags( ones(N,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], N, N );
Dx2(N,N-1) = 2/dx^2;


% solve for u0

temp = 3*mu([0;0;th0]);     % add ghost node to th
tiph = ( temp(1:end-1) + temp(2:end) ) ./ 2;
A0  = [ D; A0  ];           % add ghost node to A
% A0   = ones(N+1,1);
% tiph = ones(N+1,1);






% Dx2u = spdiags( [ A0(2:end).*tiph(2:end), -(A0(1:end-1).*tiph(1:end-1)+A0(2:end).*tiph(2:end)), A0(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], N, N );
% Dx2u(1,2) = Dx2u(1,2) + A0(1).*tiph(1) / dx^2;              % include effect from Neumann BC
% fu   = - St*( A0(1:end-1) + A0(2:end) )/ 2;
% %fu(1)  = fu(1)   - 2*P0/(3*D*dx); % include derivative (again, Neumann BC)
% fu(1) = fu(1) -2*tiph(1)*A0(1)*P0/(D*dx); 
% 
% %fu(1) =  fu(1) - 2*(1/A0(1)*(A0(3)-A0(2))/(dx)); 
% fu(end)= fu(end) - 1   * A0(end)* tiph(end)/dx^2;
% 
% u(:,1) = Dx2u\fu;

u(:,1)=1./A(:,1);


for i=2:K
    %% Solve for A at next time step (explicit, hyperbolic)
    if sum( abs(u(:,i)) > dx/dt )
        err('CFL condition broken!')
    end
    
%     tmpU = [1/D; u(:,i-1)];
    tmpU = [ u(:,i-1); 1 ];
    tmpA = [2*D-A(1,i-1); A(:,i-1);1];
    FL = tmpA(1:end-1).*tmpU;
    FR = tmpA(2:end  ).*tmpU;
    F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) ./ 2 ;
    F  = ( F(1:end-1) - F(2:end) ) ./ dx;
    A(:,i) = A(:,i-1) + dt .* F;
    A(1,i) = D;
    
    % possibly correct if A>1 somewhere (open can of worms)
        
        
    %% Solve for th at next time step (semi-implicit, parabolic)
    % Assemble matrices
    % - first order derivative
    temp = [ A(:,i); 1 ];                                
    g    = 2./( A(:,i) + [A(2:end,i);1] ) .* ( temp(2:end) - temp(1:end-1) )/dx - tmpU(2:end);  % notice u is treated explicitly
    Dx1 = spdiags( (ones(N,1).*g).*[ 1, -1 ]/(2*dx), [-1,1], N, N )';
    Dx1(N,N-1) = 0;
    % - reaction term
    Z = 2*Bi/Pe * spdiags( 1./sqrt(([A(2:end,i);1]+A(:,i))/2), 0, N, N );
    % - assemble matrix
    M = speye(N) - dt*( Dx1 + 1/Pe * Dx2 - Z );
    % - assemble rhs
     fth = (2*Bi/Pe ./sqrt(A(:,i-1))).*tha + q;   % CHANGE ME BACK TOO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%    fth =  q;
    % - solve
    th(:,i) = M\( th(:,i-1) + dt*fth );
    
    
    
    %% Solve for u at next time step (laplacian)
    temp = 3*mu([0;0;th(:,i)]);     % add ghost node to th
    tiph = ( temp(1:end-1) + temp(2:end) ) ./ 2;
    A0  = [ D; A(:,i)  ];           % add ghost node to A

    Dx2u = spdiags( [ A0(2:end).*tiph(2:end), -(A0(1:end-1).*tiph(1:end-1)+A0(2:end).*tiph(2:end)), A0(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], N, N );
    Dx2u(1,2) = Dx2u(1,2) + A0(1).*tiph(1) / dx^2;              % include effect from Neumann BC
    fu   = - St.*( A0(1:end-1) + A0(2:end) )./ 2;

    fu(1) = fu(1) -2*tiph(1)*A0(1)*P0/(D*dx);  % include derivative (again, Neumann BC)
    
    fu(end)= fu(end) - 1   * A0(end)* tiph(end)/dx^2;

    u(:,i) = Dx2u\fu;

    
    
%     Dx2u = spdiags( [ A(2:end,i).*tiph(2:end), -(A(1:end-1,i).*tiph(1:end-1)+A(2:end,i).*tiph(2:end)), A(1:end-1,i).*tiph(1:end-1) ] / dx^2, [-1,0,1], N-1, N-1 );
%     Dx2u(1,2) = Dx2u(1,2) + A0(1).*tiph(1) / dx^2;              % include effect from Neumann BC
%     fu   = St*( A(1:end-1,i) + A(2:end,i) )/ 2;
%     fu(1)  = fu(1)   - 1/D * A(1,i)* tiph(1)  /dx^2;
%     fu(end)= fu(end) - 1   * A(N,i)* tiph(end)/dx^2;
% 
%     u(1:N-1,i) = Dx2u\fu;
%     u(N,i)     = 1;
 
end


th = [ zeros(1,K); ...
       th          ];

A  = [ D*ones(1,K); ...
       A           ];

% u  = [ 1/D*ones(1,K); ...
%        u           ];
u  = [ u ; ...
       ones(1,K)   ];
   
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

