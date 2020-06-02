%function [ th, A, u, x, t ] = TimeDependentFDfullMOL( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N,uf,plt)
function yt = coupledPde(t,y)
% One vector to two vectors
global Pe Bi tha L K D uf x1 x2


%% Extract the values from the vector
A = y(1:K); 
w = y(K+1:2*K); 
th = w./A; 
phi = y(2*K+1:3*K);

% phi is valid for x \in (lambda, 1), that is , X \in (1, 0). note that
% this is in reverse order. The initial condition is given in this reverse
% order ! So in theory, I do not need to flip anything, I will just
% calculate it like this
lam = y(3*K+1); 
%FK  = y(3*K+2); 

Qvalue = 1;
Q = @(x) Qvalue*(x>x1).*(x<x2);    % heat source  
%Qth = @(x) Qvalue*(x>x1/lam).*(x<x2/lam);   % heat source
%Qphi = @(x) Qvalue*(x<((L-x1)./(L-lam))).*(x>((L-x2)./(L-lam)));    % heat source at the bottom

%% Finite volumes for A 

dX = 1/(K-1);
X = linspace(0,1,K)'; % Since we have Dirichlet boundary conditions, we don't need x=0

%A(end) = 1; 
u = usolution(A,th,lam,1); 
% This is the one that works. If you want to get it to work as before,
% revert back to this one, K
lamt = 1+(1-u(end))./(2*(1-A(end)));  % compute lamt with the edges
%lamt2 = 1+2*(1-u(end-1))./((A(end)-A(end-1)));  % compute lamt with the edges
%lamt2 = 1+(1-u(end-1))./((A(end)-A(end-1)));  % compute lamt with the edges
lamt2 = lamt; 
tmpU = [ u; uf ] -lamt.*[X; 1] ; % size K+1 x 1 

% add ghost node to A
tmpA = [D; A; 1]; % extend A by 1
S = -A.*lamt./(lam); 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;

F  = ( F(1:end-1) - F(2:end) ) / dX; % F is the rhs to dA/dt

Arhs = F./lam+S; 


 %% Solve for w at next time step with finite volumes  (semi-implicit, parabolic)
 % Assemble matrices
% Calculate fluxes for w 
phi(end) = w(end); 
% Add ghost node to w (and extend by one term)
tmpw = [0; w; phi(end)]; %

FLw = tmpw(1:end-1).*tmpU; % we use the same velocity as for A
FRw = tmpw(2:end  ).*tmpU; 
% Add ghost node to A and ghost node to th (and extend by one term)
Atmp  = [ D; A];  
thtmp = [ 0; th; phi(end)];

% Calculate flux
Fw = (FLw + FRw - abs(tmpU).*(tmpw(2:end) - tmpw(1:end-1)) ) / 2 +  ...
    (1./((dX).*Pe.*(lam))).*(- Atmp.*(thtmp(2:end)-thtmp(1:end-1))); 

Fw = (Fw(1:end-1) - Fw(2:end) )./(lam.* dX); 

% Calculate source term 
Aw = (Atmp(1:end-1) + Atmp(2:end))/2;

Sw = - (w.*lamt)./lam  - ...
    (2*Bi./Pe).*((Aw).^(1/2)).*(th-tha) + Q(X(1:end).*lam).*(Aw);

%% Solve for phi at next time step with finite volumes  (semi-implicit, parabolic)
% Calculate fluxes for phi
Xbar = linspace(1,0,K+1)';

dXbar = 1/(K-1); 
%lamt2 = 1+(1-u(end-1))./((A(end)-A(end-1)));  % compute lamt with the edges

tmpU = lamt2.*Xbar -ones(size(Xbar)); 

tmphi = [phi(1); phi; phi(end)]; % 

FLp = tmphi(1:end-1).*tmpU;
FRp = tmphi(2:end  ).*tmpU; 


Fp = (FLp + FRp + abs(tmpU).*(tmphi(1:end-1) - tmphi(2:end) ) ) / 2 - ...
    (tmphi(2:end)-tmphi(1:end-1)) ./ (dXbar.*Pe*(L-lam)); 

% Calculate flux differences
Fp = (Fp(1:end-1) - Fp(2:end) ) ./((L-lam).*dXbar) ; 

% Calculate source terms for phi
Sp =  lamt2.*phi./(L-lam)  - (2*Bi./Pe).*(phi-tha) + ...
    flip(Q(L-(Xbar(1:end-1).*(L-lam)))); 


%% This is what I had set originally

% Continuity of flux
%Fp(end) = Fw(end); 
%Sp(end) = Sw(end); 

% Fw(end) = Fp(end); 
% Sw(end) = Sp(end); 



%% Assemble RHS

wt = Fw + Sw;  
phit = Fp +Sp ;

At = Arhs; 

% - solve

% I set the last derivative of A to zero, so that A remains 1 at the end.
% (that way, when I plot A, I will get a nice continuous change from A
% almost 1 to 1. 

yt(1:K)= At;
%yt(1:K)= [At(1:end-1); 0];

%% Impose RHS

 phit(end)= wt(end);
% phit(end) = 0; 
yt(K+1:2*K) = wt;
yt(2*K+1:3*K) = phit; 
yt(3*K+1) = lamt; 
%yt(3*K+2) = FK; 
yt = yt'; 

end