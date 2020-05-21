%function [ th, A, u, x, t ] = TimeDependentFDfullMOL( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N,uf,plt)
function yt = coupledPde(t,y)
% One vector to two vectors
global Pe Bi tha L K D uf x1 x2
% A, w, and th are valid for x \in (0,lambda), that is, X \in (0, 1).
A = y(1:K); % This always imports the initial condition, so I don't really understand where to change u...
w = y(K+1:2*K); 
th = w./A; 

% lambda is constant in space ( so, x, X, Xbar give the same value)
lam = y(2*K+1:3*K);

% phi is valid for x \in (lambda, 1), that is , X \in (1, 0). note that
% this is in reverse order. The initial condition is given in this reverse
% order ! So in theory, I do not need to flip anything, I will just
% calculate it like this
phi = y(3*K+1:4*K); 
%th(end) = w(end); 

% heat source: this is a function in x, meaning that if I want to evaluate
% it at X, I need to send X*lam; similarly, if I want to evaluate it at
% Xbar, I need to send L - Xbar(L-lamb). 

Qvalue = 1;
Q = @(x) Qvalue*(x>x1).*(x<x2);    % heat source  
%Qth = @(x) Qvalue*(x>x1/lam(end)).*(x<x2/lam(end));   % heat source
%Qphi = @(x) Qvalue*(x<((L-x1)./(L-lam(end)))).*(x>((L-x2)./(L-lam(end))));    % heat source at the bottom
%% We solve everything for the TOP 
% We extract the necessary vectors 

dX = 1/K;
X = ( dX:dX:1 )'; % Since we have Dirichlet boundary conditions, we don't need x=0
% solve for u
A(end) = 1; 
u = usolution(A,th,lam(end),1); 
%tmpA_edge =  (A + [A(2:end);1])/2; % extract A at the edges 
%lamt = 1+(u(end)-u(end-1))./(1-tmpA_edge(end-1));  % compute lamt with the edges
lamt = 1+(u(end)-u(end-1))./(A(end)-A(end-1));  % compute lamt with the edges
tmpU = [ u; uf ] -lamt.*[X; 1] ; % size K+1 x 1 
%A(end) = 1; 
% add ghost node to A
tmpA = [D; A; 1]; % extend A by 1
S = -A.*lamt./(lam); 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;

F  = ( F(1:end-1) - F(2:end) ) / dX; % F is the rhs to dA/dt
%F(end) = 1; S(end) = 0 ; 
Arhs = F./lam+S; 


 %% Solve for th at next time step with finite volumes  (semi-implicit, parabolic)
 % Assemble matrices
% Calculate fluxes for w 
%tmpw = [0; w; A(end)*th(end)]; %
tmpw = [0; w; phi(end)]; %
FLw = tmpw(1:end-1).*tmpU;
FRw = tmpw(2:end  ).*tmpU; 
Atmp  = [ D; A; 1];  
thtmp = [ 0; th; phi(end)]; % edges or no edges, it is a similar solution 

Fw = (FLw + FRw + abs(tmpU).*(tmpw(1:end-1) - tmpw(2:end) ) ) / 2 +  ...
    (1./((dX).*Pe.*(lam(end)))).*(- Atmp(1:end-1).*(thtmp(2:end)-thtmp(1:end-1))); 
Fw = (Fw(1:end-1) - Fw(2:end) )./(lam(end).* dX); 
%thtmp = (thtmp(1:end-1)+thtmp(2:end))/2;   % extract th at the edges    

%Sw = - (A.*lamt.*thtmp(2:end))./lam  - ...
 %   (2*Bi./(Pe)).*sqrt(([A(2:end);1]+A)/2).*(thtmp(2:end)-tha) + Q(X*lam(end)).*(([A(2:end);1]+A)/2);
Sw = - (A.*lamt.*th)./lam  - ...
    (2*Bi./(Pe)).*sqrt(([A(2:end);1]+A)/2).*(th-tha) + Q(X*lam(end)).*(([A(2:end);1]+A)/2);

% HERE, WE WANT TO SOLVE, FROM X=1 TO X=0, 
% Calculate fluxes for phi
Xbar = linspace(1,0,K)';
dXbar = 1/K; 
%lamt = -lamt; 
tmpU = lamt.*[1; Xbar] -1; % the scaling of U is outside of the flux function. size K+1 x 1 THIS IS U, FROM X=1 TO X=0  [lamt should be constant anyway]
% Here, phi(1) corresponds to phi at x=1, X=0, and phi(end) corresponds to
% phi at x=lam, X=1. 
phi(end) = th(end); 
tmphi = [phi(1); phi; th(end)]; % In theory, this "1" shouldn't matter at all,  I think I am imposing neumann here ... by accident
% tmphi = flip(tmphi); % phi was provided flipped ... 
FLp = tmphi(1:end-1).*tmpU;
FRp = tmphi(2:end  ).*tmpU; 
%tmphi_edge = (phi + [phi(2:end); th(end)])/2;
%tmphi_edge = (tmphi(1:end-1)+tmphi(2:end))/2;   % extract th at the edges
%tmphi_edge = [phi(1); tmphi_edge; th(end)]; 
%phix = (tmphi(3:end) - tmphi(1:end-2))./(2*dXbar);
phix = derivative(phi,dXbar)'; 
phix= [0; phix]; % ????? does this make sense?  
%phix(1) = 0; 
%phix(end) = (th(end)-tmphi_edge(end-1))/dXbar; 
Fp = (FLp + FRp + abs(tmpU).*(tmphi(1:end-1) - tmphi(2:end) ) ) / 2 - (phix) ./ (Pe*(L-lam(end))); 

% Calculate flux differences
Fp = (Fp(1:end-1) - Fp(2:end) ) ./((L-lam(end)).*dXbar) ; 
  
Fp(end) = Fw(end); % If ( change this, it doesn't matter what phix end is. 

% Calculate source terms
Sp =  lamt.*phi./(L-lam(end))  - (2*Bi./Pe).*(phi-tha) + Q(L-flip(Xbar).*(L-lam(end))); 
%  plot(Qphi(Xbar))
%  hold on 
 % Calculate the rhs
wt = Fw + Sw;  
phit = Fp +Sp ;

% - assemble rhs
At = Arhs; 
% - solve
%yt(1:K)= At;
yt(1:K)= [At(1:end-1); 0];

% Don't set these guys to zero as this creates a maximum at that point,
% given that the derivative wrt x is also zero (or close to it) 
wt(end) = 0;
yt(K+1:2*K) = wt;
phit(end) = 0; 
yt(2*K+1:3*K) = lamt; 
yt(3*K+1:4*K) = phit; 
yt = yt'; 

end