%function [ th, A, u, x, t ] = TimeDependentFDfullMOL( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N,uf,plt)
function yt = coupledPdeNoFBwithFV(t,y)
% One vector to two vectors
global Pe Bi tha L K D uf 
% for i=1:n
%     A(i) = y(i); 
%     th(i) = y(i+n);
% end
dx = L/K;
x = ( dx:dx:L )'; % Since we have Dirichlet boundary conditions, we don't need x=0

A = y(1:K); % This always imports the initial condition, so I don't really understand where to change u...
w = y(1+K:2*K);

th = w./A; % I think that this is okay 
% viscosity
%mu = @(th) exp(-gamma*th); % Perhaps discretise this as well, but put a pin in that

% heat source
x1 = 5/7;
x2 = 6.5/7;
Qvalue = 1;
Q = Qvalue*(x>x1).*(x<x2);    % heat source  
%q = zeros(N,1);


u = usolutionNoFB(A,th); 

tmpU = [ u; uf ];

% add ghost node to A

tmpA = [D; A; 1]; % extend A by 1
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;
F  = ( F(1:end-1) - F(2:end) ) / dx; % F is the rhs to dA/dt

At = F; 
  %% Solve for th at next time step with finite volumes  (semi-implicit, parabolic)
 % Assemble matrices
% opt = 1;  
%w = A.*th;
%if opt == 1 
    tmpw = [0; w; w(end)]; % In theory, this "1" shouldn't matter at all,  I think I am imposing neumann here ... by accident
    FL = tmpw(1:end-1).*tmpU;
    FR = tmpw(2:end  ).*tmpU; 
    F = (FL + FR + abs(tmpU).*(tmpw(1:end-1) - tmpw(2:end) ) ) / 2 ; 
     
    Atmp  = [ D; A ];  

    thtmp = [ 0; th; th(end)]; % edges or no edges, it is a similar solution 
    thtmp = (thtmp(1:end-1)+thtmp(2:end))/2;   % extract th at the edges
    
    F = (F(1:end-1) - F(2:end) ) / dx +(1./((dx^2).*Pe)).*(Atmp(2:end).*([thtmp(3:end); thtmp(end)]- thtmp(2:end))); ; 
   
   
   % S = (1./((dx^2).*Pe)).*(Atmp(2:end).*([thtmp(3:end); thtmp(end)]- thtmp(2:end))  - Atmp(1:end-1).*(thtmp(2:end)-thtmp(1:end-1))) - ...
   %     (2*Bi./(Pe)).*sqrt(([A(2:end);1]+A)/2).*(thtmp(2:end)-tha) + Q.*(([A(2:end);1]+A)/2); 
     S = (1./((dx^2).*Pe)).*(- Atmp(1:end-1).*(thtmp(2:end)-thtmp(1:end-1))) - ...
        (2*Bi./(Pe)).*sqrt(([A(2:end);1]+A)/2).*(thtmp(2:end)-tha) + Q.*(([A(2:end);1]+A)/2); 

% mat = spdiags( [Atmp(2:end), -(Atmp(1:end-1) + Atmp(2:end)), Atmp(1:end-1)] / dx^2 , [-1 0 1], K,K ); 
% vec = mat*thtmp(2:end); 
% 
% S2 = (1/Pe).*vec -(2*Bi./(Pe)).*sqrt(([A(2:end);1]+A)/2).*(thtmp(2:end)-tha) + Q.*(([A(2:end);1]+A)/2);

wt = F + S; 

%else
    % this part below DOES NOT WORK, but should be the framework in case I
    % need to incorporate stuff into the flux
    % BUT IT DOES NOT WORK! IN FACT, IF YOU DO NOT NEED THE DETAILS BELOW
    % FOR THE FREE BOUNDARY, DO.NOT.BOTHER. JUST DELETE.
%     tmpw = [0; w; w(end)]; % In theory, this "1" shouldn't matter at all
%     FL = tmpw(1:end-1).*tmpU;  % add diffusive term
%     FR = tmpw(2:end  ).*tmpU; % add diffusive term; 
%     F = (FL + FR + abs(tmpU).*(tmpw(1:end-1) - tmpw(2:end) ) ) / 2 ; 
%     F = (F(1:end-1) - F(2:end) ) / dx - (1/Pe).*(tmpw(2:end-1)-tmpw(1:end-2)./dx) ; 
%     %thtmp = (th + [ th(2:end); th(end)])/2; 
%     thtmp = [ 0; th]; % edges or no edges, it is a similar solution 
%     %thtmp = (thtmp(1:end-1)+thtmp(2:end))/2;   % extract th at the edges
%     %thtmp = [thtmp; thtmp(end)]; % extend edge by one 
%     Anodes = ([A(2:end);1]+A)/2;
%     Anodes = [D; Anodes; 1]; 
%     
%     S = (1/Pe).*(thtmp(1:end-1).*(Anodes(2:end-1)-Anodes(1:end-2)) - thtmp(2:end).*(Anodes(3:end)-Anodes(2:end-1)))./(dx^2) - ...
%     (2*Bi./(Pe)).*(1./sqrt(Anodes(2:end-1)).*(w-tha) )+ Q.*Anodes(2:end-1); 
%     wt = F + S; 

%end 


% - solve
yt(1:K)= At;
yt(K+1:2*K) = wt; 
yt = yt'; 

end

