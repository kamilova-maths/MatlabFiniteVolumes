%% CODE THAT DETERMINES A AND LAMBDA(t), NO TEMPERATURE ADDED
% the viscosity is constant and chosen such that A0 reaches 1 at around x=1
% Each variable is N,K, where N is discretisation in t and K is
% discretisation in x 

function yt = coupledPdeNoTemp(t,y)
% One vector to two vectors
global L K gamma P0 St D uf 
% for i=1:n
%     A(i) = y(i); 
%     th(i) = y(i+n);
% end
dx = L/K;
x = ( dx:dx:L )'; % Since we have Dirichlet boundary conditions, we don't need x=0

A = y(1:K); % This always imports the initial condition, so I don't really understand where to change u...
% th = y(1+K:2*K);
% u = y(2*K+1:3*K); 
lam = y(1+K:2*K); 


theta = zeros(size(A)); 
%u = y(2*K+1:3*K); 


% viscosity
%mu = @(th) exp(-gamma*th); % Perhaps discretise this as well, but put a pin in that

% heat source
x1 = 5/7;
x2 = 6.5/7;
Qvalue = 1;
Q = Qvalue*(x>x1).*(x<x2);    % heat source  
%q = zeros(N,1);

% solve for u
u = usolution(A,theta,lam(end)); 


% Original that works sort of okay 
tmpA =  (A + [A(2:end);1])/2; % extract A at the edges 
lamt = 1+(u(end)-u(end-1))./(tmpA(end)-tmpA(end-1));  % compute lamt with the edges

% % Other iteration which does not work 
%lamt = 1+(u(end)-u(end-1))./(A(end)-A(end-1));  % BIG NO NO 
%lamt = 1+2*(u(end)-u(end-1))./(A(end)-A(end-1));  % ANOTHER BIG NO NO 


%  Atmp  = [ D; A ];   
%  Atmp  =  (Atmp(1:end-1) + Atmp(2:end))/2; % extract A at the edges 
%  lamt  = 1+(u(end)-u(end-1))./(Atmp(end)-Atmp(end-1));  % compute lamt with the edges

%tmpU = [ u; uf ];
%tmpU = [ u; uf ] -lamt.*[x; 1] ; % size K+1 x 1 
tmpU = [ u; uf ] -lamt.*[0; x] ; % size K+1 x 1 

% Source term


% add ghost node to A
tmpA = [D; A; 1]; % extend A by 1

%S  = -(A + [A(2:end);1]).*lamt./(2*lam); 
S = -A.*lamt./(lam(end));  % this (or the above expression) makes no difference, 
% stop changing it back and forth 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;
F  = ( F(1:end-1) - F(2:end) ) / dx; % F is the rhs to dA/dt
Arhs = F./lam(end)+S; 


At = Arhs; 
% - solve
yt(1:K)= At;

yt(K+1:2*K) = lamt ; 

yt = yt'; 

end

