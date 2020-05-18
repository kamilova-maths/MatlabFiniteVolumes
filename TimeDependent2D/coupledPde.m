%function [ th, A, u, x, t ] = TimeDependentFDfullMOL( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N,uf,plt)
function yt = coupledPde(t,y)
% One vector to two vectors
global Pe Bi tha L K D uf 
% for i=1:n
%     A(i) = y(i); 
%     th(i) = y(i+n);
% end
dx = L/K;
x = ( dx:dx:L )'; % Since we have Dirichlet boundary conditions, we don't need x=0

A = y(1:K); % This always imports the initial condition, so I don't really understand where to change u...
w = y(K+1:2*K); 
th = w./A; 
lam = y(2*K+1:3*K);

phi = y(3*K+1:4*K); 

% heat source
x1 = 5/7;
x2 = 6.5/7;
Qvalue = 1;
Q = @(x) Qvalue*(x>x1).*(x<x2);    % heat source  
Qth = Qvalue*(x>x1/lam(end)).*(x<x2/lam(end));   % heat source
Qphi = Qvalue*(x<((1-x1)./(1-lam(end)))).*(x>((1-x2)./(1-lam(end))));    % heat source at the bottom
%% We solve everything for the TOP 
% We extract the necessary vectors 
% solve for u
u = usolution(A,th,lam(end)); 
tmpA =  (A + [A(2:end);1])/2; % extract A at the edges 
lamt = 1+(u(end)-u(end-1))./(tmpA(end)-tmpA(end-1));  % compute lamt with the edges
tmpU = [ u; uf ] -lamt.*[x; 1] ; % size K+1 x 1 
% add ghost node to A
tmpA = [D; A; 1]; % extend A by 1
S = -A.*lamt./(lam); 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;

F  = ( F(1:end-1) - F(2:end) ) / dx; % F is the rhs to dA/dt
Arhs = F./lam+S; 


 %% Solve for th at next time step with finite volumes  (semi-implicit, parabolic)
 % Assemble matrices
% Calculate fluxes for w 
%tmpw = [0; w; A(end)*th(end)]; %
tmpw = [0; w; th(end)]; %
FLw = tmpw(1:end-1).*tmpU;
FRw = tmpw(2:end  ).*tmpU; 
Atmp  = [ D; A; 1];  
thtmp = [ 0; th; th(end)]; % edges or no edges, it is a similar solution 
Fw = (FLw + FRw + abs(tmpU).*(tmpw(1:end-1) - tmpw(2:end) ) ) / 2 +  ...
    (1./((dx).*Pe.*(lam(end)))).*(- Atmp(1:end-1).*(thtmp(2:end)-thtmp(1:end-1))); 

% Calculate fluxes for phi
tmpU = lamt.*[x; 1] -1 ; % size K+1 x 1 
tmpU = flip(tmpU); 
tmphi = [phi(1); phi; th(end)]; % In theory, this "1" shouldn't matter at all,  I think I am imposing neumann here ... by accident
tmphi = flip(tmphi); 
FLp = tmphi(1:end-1).*tmpU;
FRp = tmphi(2:end  ).*tmpU; 
tmphi_edge = (tmphi(1:end-1)+tmphi(2:end))/2;   % extract th at the edges

phix = derivative(tmphi_edge,dx)'; 
phix(end) = 0; 
Fp = (FLp + FRp + abs(tmpU).*(tmphi(1:end-1) - tmphi(2:end) ) ) / 2 - (phix) ./ (Pe*(1-lam(end))); 

% Calculate flux differences

Fw = (Fw(1:end-1) - Fw(2:end) )./(lam(end).* dx); 
Fp = (Fp(1:end-1) - Fp(2:end) ) ./((1-lam).*dx) ; 
  
%Fp(1) = Fw(end); 

% Calculate source terms
thtmp = (thtmp(1:end-1)+thtmp(2:end))/2;   % extract th at the edges    
Sw = - (A.*lamt.*thtmp(2:end))./lam  - ...
    (2*Bi./(Pe)).*sqrt(([A(2:end);1]+A)/2).*(thtmp(2:end)-tha) + Qth.*(([A(2:end);1]+A)/2);

Sp =  (lamt.*tmphi_edge(2:end))./(1-lam)  - (2*Bi./(Pe)).*(tmphi_edge(2:end)-tha) + flip(Qphi); 
 
 % Calculate the rhs
wt = Fw + Sw;  
phit = Fp +Sp ;

%flip back phi
phit = flip(phit); 


% - assemble rhs
At = Arhs; 

% - solve
yt(1:K)= At;
yt(K+1:2*K) = wt;  
yt(2*K+1:3*K) = lamt; 
yt(3*K+1:4*K) = phit; 
yt = yt'; 

end

