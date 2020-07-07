%% CODE THAT DETERMINES A AND LAMBDA(t), NO TEMPERATURE ADDED
% the viscosity is constant and chosen such that A0 reaches 1 at around x=1
% Each variable is N,K, where N is discretisation in t and K is
% discretisation in x 

function yt = coupledPdeNoTemp(t,y)
% One vector to two vectors
global L K P0 D uf 



A = y(1:K); % This always imports the initial condition, so I don't really understand where to change u...
% th = y(1+K:2*K);
% u = y(2*K+1:3*K); 
lam = y(1+K); 


theta = zeros(size(A)); 


X = linspace(0,1,K+1)'; % Since we have Dirichlet boundary conditions, we don't need x=0
dX = 1/K; 
dXbar = 1/K;

u = usolution(A,theta,lam,1,P0); 
dudx = (1/2).*(uf-u(end))/(lam*dX);
dAdx = 2*(1-A(end))/(lam*dX + (L-lam)*dXbar);
lamt = 1+(dudx/dAdx);
%lamt = 1+(uf-u(end))/(2*(1-A(end)));  % compute lamt with the edges

tmpU = [ u; uf ] -lamt.*X ; % size K+1 x 1 

%tmpU = [ 1/D; u ] -lamt.*X ; % size K+1 x 1 


% add ghost node to A
%tmpA = [D; A; 1]; % extend A by 1
tmpA = [2*D-A(1); A; 1]; 
%S  = -(A + [A(2:end);1]).*lamt./(2*lam); 
S = -A.*lamt./(lam);  % this (or the above expression) makes no difference, 
% stop changing it back and forth 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
%exc = tmpU(end).*tmpA(end)+lamt
%S(end) = S(end)+lamt; 
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;
F  = ( F(1:end-1) - F(2:end) ) / (dX*lam); % F is the rhs to dA/dt
Arhs = F +S; 
Arhs(end) = Arhs(end)+lamt; 

At = Arhs; 
% - solve
yt(1:K)= At;

yt(K+1) = lamt ; 

yt = yt'; 

end

