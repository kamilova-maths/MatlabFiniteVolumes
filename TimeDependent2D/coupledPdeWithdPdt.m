function yt = coupledPdeWithdPdt(t,y)
% DATE:  2020
%
% DESCR:    yt = coupledPdeWithPdt(t,y)
%           Code that encapsulates the partial differential equations for 
%           the time dependent version of the flow and temperature problem, 
%           where we additionally solve for dPdt. There are no events
%           associated with this function.
%           We use the method of lines here, so we construct the matrices
%           for the x discretisation here, whereas the t discretisation is
%           being done automatically in the main code by ode15s.
%           Additionally, ode15s.m invokes this code. 
%
% INPUT: 
%           t:  Particular t value using by the internal MATLAB routine
%           ode15s.m
%           y:  Long vector that includes all the variables of interest. It
%           is of the same length as y0, the initial condition vector,
%           provided in the main code that invokes ode15s.m
% OUTPUT:   yt: Vector of the time derivatives for each variable of
%           interest. 
% ADDITIONAL COMMENTS: 
%           This code uses should be run for a certain amount of time or 
%           until some condition is satisfied so that we obtain some kind
%           of balance in the dPdt equation. Otherwise, P will just
%           decrease over time and become negative. This code should be
%           used in control systems designed to model the lowering of paste
%           cylinders into electrode casings.
%
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           bvpinit, bvp4c: to solve the ODEs
%           usolution: Routine that solves the u problem, which is
%           independent of time, so it does not need to be in the time
%           vector. We compute the form of u at each timestep t to use it
%           in the temperature and area equations. 
%           EventFunction: Checks if P is smaller than some set tolerance,
%           and interrupts the computation if that's the case.

% One vector to two vectors
global Pe Bi tha L K D uf x1 x2 Q St


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
u = usolution( A, th, lam, 1, P);

Pt = D*St.*(- u(1));

lamt = 1+(1-u(end))/(2*(1-A(end)));  % compute lamt with the edges


tmpU = [u; uf ] -lamt.*X ; % size K+1 x 1 


tmpA = [2*D-A(1); A; 1]; % extend A by 1

 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;



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