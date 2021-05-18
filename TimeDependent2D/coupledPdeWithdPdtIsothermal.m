function yt = coupledPdeWithdPdtIsothermal(t,y)
% DATE:  2020
%
% DESCR:    yt = coupledPdeWithdPdtIsothermal(t,y)
%           Code that encapsulates the partial differential equations for 
%           the time dependent version of the flow problem in the isothermal
%           case . We impose a time dependent equation for P as well, and
%           we require another code where this is done inside a loop, so
%           that the code is restarted. Otherwise, the solution will never
%           reach a steady state.
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
%           This code uses P0t, a time dependent version of P0. If we want 
%           to use a constant P0t, the best way is to set the time
%           dependent function to be solely defined by P0.  
%
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           bvpinit, bvp4c: to solve the ODEs
%           usolution: Routine that solves the u problem, which is
%           independent of time, so it does not need to be in the time
%           vector. We compute the form of u at each timestep t to use it
%           in the temperature and area equations. 



global K D uf St

%% Extract the values from the vector
Alamt = y(1:K); 

P = y(K+1); 
lam = y(K+2);

A = Alamt./lam;

%% discretisation parameters - these correspond to the interfaces. 
% - left part of domain [0,\lambda]
dX = 1/K;
X = linspace(0,1,K+1)'; % Since we have Dirichlet boundary conditions, we don't need x=0

%% Finite volumes for A 
u = usolution( A, zeros(size(A)), lam, 1, P);

Pt = D*St.*(- u(1));

lamt = 1+(uf -u(end))/(2*(1-A(end)));  % compute lamt with the edges

tmpU = [u; uf ] -lamt.*X ; % size K+1 x 1 

tmpA = [2*D-A(1); A; 1]; % extend A by 1
 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;

F  = ( F(1:end-1) - F(2:end) ) / (dX); % F is the rhs to dA/dt

Alamtrhs = F; 

 

%% Assemble RHS

Alamt = Alamtrhs; 

%% Impose RHS
yt(1:K)= Alamt;
yt(K+1) = Pt; 
yt(K+2) = lamt; 
yt = yt'; 

end