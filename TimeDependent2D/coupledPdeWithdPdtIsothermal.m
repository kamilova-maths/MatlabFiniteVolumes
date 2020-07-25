function yt = coupledPdeWithdPdtIsothermal(t,y)
% One vector to two vectors
global K D uf St

%% Extract the values from the vector
Alamt = y(1:K); 

P = y(K+1); 
lam = y(K+2);

A = Alamt./lam;

%% discretisation parameters - these correspond to the interfaces. 
% - left part of domain [0,\lambda]
dX = 1/K;
X = linspace(0,1,K+1)'; % Since we have Dirichlet boundary conditions, we don't need x=0

%% Finite volumes for A 
u = usolution( A, zeros(size(A)), lam, 1, P);

Pt = D*St.*(- u(1));

lamt = 1+(uf -u(end))/(2*(1-A(end)));  % compute lamt with the edges

tmpU = [u; uf ] -lamt.*X ; % size K+1 x 1 

tmpA = [2*D-A(1); A; 1]; % extend A by 1
 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;

% This works but doesn't change the slow changing drift. I think my other
% control system is better anyway 

%F(end) = 1 - lamt;

F  = ( F(1:end-1) - F(2:end) ) / (dX); % F is the rhs to dA/dt

Alamtrhs = F; 

 

%% Assemble RHS

Alamt = Alamtrhs; 

%% Impose RHS
yt(1:K)= Alamt;
yt(K+1) = Pt; 
yt(K+2) = lamt; 
yt = yt'; 

end