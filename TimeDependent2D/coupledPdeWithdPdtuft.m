%function [ th, A, u, x, t ] = TimeDependentFDfullMOL( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N,uf,plt)
function yt = coupledPdeWithdPdtuft(t,y)
% One vector to two vectors
global Pe Bi tha L K D x1 x2 Q St uft


%% Extract the values from the vector
Alamt = y(1:K); 

w = y(K+1:2*K); 

phi = y(2*K+1:3*K);
P = y(3*K+1); 
lam = y(3*K+2);

A = Alamt./lam;
th = w./A; 
Qfun = @(x) Q*(x>x1).*(x<x2);    % heat source  

%% discretisation parameters - these correspond to the interfaces. 
% - left part of domain [0,\lambda]
dX = 1/K;
X = linspace(0,1,K+1)'; % Since we have Dirichlet boundary conditions, we don't need x=0

% - right part of domain [\lambda, L]
dXbar = 1/K;
Xbar = linspace(0,1,K+1)';

%% Finite volumes for A 
u = usolutionuft( A, th, lam, 1, P, uft(t));

%lamt = 1+(dudx/dAdx);

Pt = D*St.*(- u(1));

%Pt = D*St.*(uval- u(1));

%Pt = D*St.*(uin(t)- u(1))
%lamt = 1+(1./Aint(end) -u(end))/(2*(1-A(end)));  % compute lamt with the edges
lamt = uft(t) +(uft(t) - u(end))/(2*(1-A(end))) ; 
%lamt = 1+(1-u(end))/(2*(1-A(end)));  % compute lamt with the edges


%tmpU = [u; uf ] -lamt.*X ; % size K+1 x 1 
tmpU = [u; uft(t) ] -lamt.*X ; % size K+1 x 1 %tmpU = [u; 1./Aint(end) ] -lamt.*X ; % size K+1 x 1 
% add ghost node to A

tmpA = [2*D-A(1); A; 1]; % extend A by 1
%tmpA = [D; A; 1];
 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;

% This works but doesn't change the slow changing drift. I think my other
% control system is better anyway 
%F(end) = 1 - lamt;

F  = ( F(1:end-1) - F(2:end) ) / (dX); % F is the rhs to dA/dt

Alamtrhs = F; 

 %% Solve for w at next time step with finite volumes  (semi-implicit, parabolic)
 % Assemble matrices
% Calculate fluxes for w 

% Add ghost node to w (and extend by one term)
tmpw = [0; w; phi(1)]; % NOTE: I took the temp value from phi (continuity)

FLw = tmpw(1:end-1).*tmpU; % we use the same velocity as for A
FRw = tmpw(2:end  ).*tmpU; 
% Add ghost node to A and ghost node to th (and extend by one term)
Atmp  = ([ 2*D-A(1); A] + [A;1] ) / 2;  
%Atmp = ([D; A] + [A;1]) / 2; 
thtmp = [ 0; th; phi(1)];% NOTE: I took the temp value from phi (continuity)

% Calculate flux
Fw = (FLw + FRw - abs(tmpU).*(tmpw(2:end) - tmpw(1:end-1)) ) / 2 -  ...
    (1./((dX).*Pe.*(lam))).*(Atmp.*(thtmp(2:end)-thtmp(1:end-1))); 
% for the last interface, the cell sizes vary: FD to recover gradient must account for this
Fw(end) = (FLw(end) + FRw(end) - abs(tmpU(end)).*(tmpw(end) - tmpw(end-1)) ) / 2 -  ...
    2./((dX*lam + dXbar*(L-lam))*Pe).*(thtmp(end)-thtmp(end-1)); % A is 1 at that interface

Fw = (Fw(1:end-1) - Fw(2:end) )./(lam.* dX); 

% Calculate source term 
Sw = - (w.*lamt)./lam  - ...
    (2*Bi./Pe).*((A).^(1/2)).*(th-tha) + Qfun(((dX/2):dX:(1-dX/2))'.*lam).*(A);

%% Solve for phi at next time step with finite volumes  (semi-implicit, parabolic)

tmpU = (Xbar - ones(size(Xbar))).*lamt + ones(size(Xbar));  
tmphi = [th(end); phi; phi(end)];% NOTE: similarly, I took the temp value from th (continuity)


FLp = tmphi(1:end-1).*tmpU;
FRp = tmphi(2:end  ).*tmpU; 

tmphi = [th(end); phi; phi(end)];% NOTE: similarly, I took the temp value from th (continuity)


Fp = (FLp + FRp - abs(tmpU).*(tmphi(2:end) - tmphi(1:end-1) ) ) / 2 - ...
    (tmphi(2:end)-tmphi(1:end-1)) ./ (dXbar.*Pe*(L-lam)); 
Fp(1) = (FLp(1) + FRp(1) - abs(tmpU(1)).*(tmphi(2) - tmphi(1) ) ) / 2 - ...
    2*(tmphi(2)-tmphi(1)) ./ ((dXbar*(L-lam) + dX*lam)*Pe); 
% Fp(1) = uga;	% nah, this won't work: it's only the laplacian part of the
% flux that's influenced by the different dX

% Calculate flux differences
Fp = (Fp(1:end-1) - Fp(2:end) ) ./((L-lam).*dXbar) ; 

% Calculate source terms for phi
Sp =  lamt.*phi./(L-lam)  - (2*Bi./Pe).*(phi-tha) + ...
(Qfun(lam+(((dXbar/2):dXbar:(1-dXbar/2))'.*(L-lam)))); 

%% Assemble RHS
wt = Fw + Sw;  
phit = Fp +Sp;
Alamt = Alamtrhs; 

% - solve
yt(1:K)= Alamt;

%% Impose RHS
yt(K+1:2*K) = wt;
yt(2*K+1:3*K) = phit; 
yt(3*K+1) = Pt; 
yt(3*K+2) = lamt; 
yt = yt'; 

end