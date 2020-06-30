%% CODE THAT DETERMINES A AND LAMBDA(t), NO TEMPERATURE ADDED
% the viscosity is constant and chosen such that A0 reaches 1 at around x=1
% Each variable is N,K, where N is discretisation in t and K is
% discretisation in x 

function yt = coupledPdeNoTempdPdt(t,y)
% One vector to two vectors
global K D uf St uin

A = y(1:K); % This always imports the initial condition, so I don't really understand where to change u...
% th = y(1+K:2*K);
% u = y(2*K+1:3*K); 
lam = y(1+K); 
P = y(2+K); 

theta = zeros(size(A)); 


X = linspace(0,1,K+1)'; % Since we have Dirichlet boundary conditions, we don't need x=0
dX = 1/K; 

u = usolution(A,theta,lam,1,P); 

% I can either calculate what u(0,t) should be, or just evaluate it at the
% first timestep. Currently these give me different things ... so probably
% trust u(1) more... 

Pt = D*St.*(uin(t)- u(1));
%u0t = (lam.*dX/(A(1)-2*D+A(1))).*P/3; % There is a real and distinct possibility that you are the reason for all my troubles
%Pt = D*St.*(uin(t)- u0t);

% I THINK THERE IS A MISTAKE HERE
lamt = 1+(uf-u(end))/(2*(1-A(end)));  % compute lamt with the edges
tmpU = [ u; uf ] -lamt.*X ; % size K+1 x 1 


% add ghost node to A
%tmpA = [D; A; 1]; % extend A by 1
tmpA = [2*D-A(1); A; 1]; 

S = -A.*lamt./(lam);  % this (or the above expression) makes no difference, 
% stop changing it back and forth 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;
F  = ( F(1:end-1) - F(2:end) ) / (dX*lam); % F is the rhs to dA/dt
Arhs = F./lam +S; 


At = Arhs; 
% - solve
yt(1:K)= At;

yt(K+1) = lamt; 

yt(K+2) = Pt; 

yt = yt'; 

end

