%function [ th, A, u, x, t ] = TimeDependentFDfullMOL( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N,uf,plt)
function yt = coupledPdeFede(t,y)
% One vector to two vectors
global Pe Bi tha L K D uf x1 x2


%% Extract the values from the vector
A = y(1:K); 
w = y(K+1:2*K); 
th = w./A; 
phi = y(2*K+1:3*K-1);

% phi is valid for x \in (lambda, 1), that is , X \in (1, 0). note that
% this is in reverse order. The initial condition is given in this reverse
% order ! So in theory, I do not need to flip anything, I will just
% calculate it like this
lam = y(3*K); 

Qvalue = 1;
Q = @(x) Qvalue*(x>x1).*(x<x2);    % heat source  
%Qth = @(x) Qvalue*(x>x1/lam).*(x<x2/lam);   % heat source
%Qphi = @(x) Qvalue*(x<((L-x1)./(L-lam))).*(x>((L-x2)./(L-lam)));    % heat source at the bottom

%% Finite volumes for A 

dX = 1/K; % NOTE: I assume we have K cells per subdomain
X = linspace(0,1,K+1)'; % Since we have Dirichlet boundary conditions, we don't need x=0

A(end) = 1; 
u = usolution(A,th,lam,1); 

lamt = 1+(1-u(end))./(2*(A(end)-A(end-1)));  % compute lamt with the edges
tmpU = [ u; uf ] -lamt.*X ; % size K+1 x 1 

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

% Add ghost node to w (and extend by one term)
tmpw = [0; w; phi(1)]; % NOTE: I took the temp value from phi (continuity)

FLw = tmpw(1:end-1).*tmpU; % we use the same velocity as for A
FRw = tmpw(2:end  ).*tmpU; 
% Add ghost node to A and ghost node to th (and extend by one term)
Atmp  = [ D; A];  
thtmp = [ 0; th; phi(1)];% NOTE: I took the temp value from phi (continuity)

% Calculate flux
Fw = (FLw + FRw - abs(tmpU).*(tmpw(2:end) - tmpw(1:end-1)) ) / 2 +  ...
    (1./((dX).*Pe.*(lam))).*(- Atmp.*(thtmp(2:end)-thtmp(1:end-1))); 

Fw = (Fw(1:end-1) - Fw(2:end) )./(lam.* dX); 

% Calculate source term 
Sw = - (w.*lamt)./lam  - ...
    (2*Bi./Pe).*((A).^(1/2)).*(th-tha) + Q(X(1:end-1).*lam).*(A);

%% Solve for phi at next time step with finite volumes  (semi-implicit, parabolic)
% Calculate fluxes for phi
dXbar = 1/(K-1);						 % NOTE: the first cell of phi is actually the last of theta, so exclude that
Xbar = linspace(dXbar,1,K)'; % similarly here: we start considering unknowns for phi from the following cell 

tmpU = lamt.*Xbar -ones(size(Xbar)); 
tmpU = -tmpU; 
tmphi = [th(end); phi; phi(end)];% NOTE: similarly, I took the temp value from th (continuity)

% NOTE: it's important (for conservation of stuff) that the flux at the
% *last* interface of th mathes that of the *first* interface of phi. This
% seems to work anyway (so probably the way the flux is computed is the
% same, at the end of the day), but if you want to be extra sure, you
% probably want to set something like FLp(1) = FRw(end), or so
FLp = tmphi(1:end-1).*tmpU;
FRp = tmphi(2:end  ).*tmpU; 

Fp = (FLp + FRp + abs(tmpU).*(tmphi(1:end-1) - tmphi(2:end) ) ) / 2 - ...
    (tmphi(2:end)-tmphi(1:end-1)) ./ (dXbar.*Pe*(L-lam)); 

% Calculate flux differences
Fp = (Fp(1:end-1) - Fp(2:end) ) ./((L-lam).*dXbar) ; 

% Calculate source terms for phi
Sp =  lamt.*phi./(L-lam)  - (2*Bi./Pe).*(phi-tha) + ...
    Q(L-(Xbar(2:end).*(L-lam))); 


%% This is what I had set originally

% Continuity of flux
%Fp(end) = Fw(end); 

%% Alternative continuity of flux
% FK = Fw(end)-Fp(end); 

Fp(1) = Fw(end);
Sp(1) = Sw(end); 

%% Assemble RHS

wt = Fw + Sw;  
phit = Fp +Sp ;

At = Arhs; 

% - solve

% I set the last derivative of A to zero, so that A remains 1 at the end.
% (that way, when I plot A, I will get a nice continuous change from A
% almost 1 to 1. 

%yt(1:K)= At;
yt(1:K)= [At(1:end-1); 0];

%% Impose RHS
yt(K+1:2*K) = wt;
yt(2*K+1:3*K-1) = phit; 
yt(3*K) = lamt; 
% yt(3*K+2) = FK; 
yt = yt'; 

end