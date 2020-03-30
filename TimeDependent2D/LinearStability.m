% Define sizes/parameters
InitialiseforFD; 

% Find steady state

% Steady states are th2, A2, u2,

th0 = th2(:,1);
A0 = A2(:,1);
u0 = u2(:,1);
N = length(th0); 

dx = L/(N-1);
Qvalue = 1;
x=linspace(0,L,N)';
Q = Qvalue*(x>x1).*(x<x2);    % heat source  

dudx = derivative(u0,dx)';
dAdx = derivative(A0,dx)';
dthdx = derivative(th0,dx)';
dudx2 = derivative(dudx,dx)';
dthdx2 = derivative(dthdx,dx)'; 

% Building the Laplacian

Dx2 = spdiags( ones(N,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], N, N );

Dx1 = sparse(1:N,1:N,1/dx*ones(N,1),N,N) ...
    + sparse(2:N,1:N-1,-1/dx*ones(N-1,1),N,N); % forward difference x derivative, homogeneous Dirichlet at x = 0

% Building the matrix 


% prefac 1 - exponential term for theta
pre1 = 3.*exp(-gamma.*th0); 
pre2 = 1./A0; 
s11 = spdiags(-dudx, 0,N,N)-u0.*Dx1;
s12 = -A0.*Dx1-spdiags(dAdx,0,N,N); 
s13 = sparse(N,N);
s21 = (dudx.*pre1).*Dx1 - spdiags((gamma.*dthdx.*dudx.*pre1), 0,N,N) + spdiags(dudx2.*pre1,0,N,N) + spdiags(St.*ones(N,1),0,N,N);
s22 = (dAdx.*pre1-gamma.*dthdx.*A0.*pre1).*Dx1 +(A0.*pre1).*Dx2; 
s23 = spdiags((gamma^2.*A0.*dudx.*dthdx.*pre1),0,N,N)-(gamma.*pre1).*Dx1.*(A0.*dudx)-(gamma.*A0.*dudx.*pre1).*Dx1; 
s31 = (1/Pe).*(dthdx.*pre2).*Dx1+spdiags((pre2.*((1/Pe).*dthdx2-dthdx.*u0 - Bi./(Pe.*sqrt(A0)).*(th0-tha)+Q')),0,N,N); 
s32 = spdiags((-dthdx.*A0.*pre2),0,N,N);
s33 = A0./Pe.*pre2.*Dx2+((1/Pe).*dAdx.*pre2-A0.*u0.*pre2).*Dx1-Dx1.*(A0.*u0.*pre2)-spdiags((2*(Bi/Pe).*sqrt(A0).*pre2),0,N,N); 

%assemble the big sparse matrix
SS = [s11,s12,s13;
      s21,s22,s23;
      s31,s32,s33];
 
  % assemble the other big sparse matrix (for the generalised eigenvalue
  % problem)
  
P1 = [speye(N,N),sparse(N,N),sparse(N,N);
      sparse(N,N),sparse(N,N),sparse(N,N);
      sparse(N,N),sparse(N,N),speye(N,N)];
  
  figure(2)
  % plot the eigenvalues
  
set(0,'DefaultAxesFontSize',12,'DefaultTextInterpreter','latex');
plot(eigs(SS,P1,21,100),'o', 'MarkerFaceColor', 'b')
set(gca,'TickLabelInterpreter','latex','fontsize',15)
xlabel('Re($\sigma$)')
ylabel('Im($\sigma$)')
