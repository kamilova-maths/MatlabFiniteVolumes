function yt = coupledPdeIsothermal(t,y)
% DATE:  2020

% DESCR:    Code that determines A and lambda(t), no temperature dependence
%           the viscosity is constant and chosen such that A0 reaches 1 at around x=1
%           Each variable is N,K, where N is discretisation in t and K is
%           discretisation in x 
%           We use the method of lines here, so we construct the matrices
%           for the x discretisation here, whereas the t discretisation is
%           being done automatically in the main code by ode15s.
%           Additionally, ode15s.m invokes this code. 

% INPUT: 
%           t:  Particular t value using by the internal MATLAB routine
%           ode15s.m
%           y:  Long vector that includes all the variables of interest. It
%           is of the same length as y0, the initial condition vector,
%           provided in the main code that invokes ode15s.m
% OUTPUT:   yt: Vector of the time derivatives for each variable of
%           interest. 
% ADDITIONAL COMMENTS: 
%           This code uses P0, a constant applied load.  

% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           bvpinit, bvp4c: to solve the ODEs
%           usolution: Routine that solves the u problem, which is
%           independent of time, so it does not need to be in the time
%           vector. We compute the form of u at each timestep t to use it
%           in the temperature and area equations. 
global K P0t D uf 


Alam = y(1:K); % Alam is the multiplication of A times lambdat
lam  = y(1+K); 
A    = Alam./lam; 


theta = zeros(size(A)); 


X = linspace(0,1,K+1)'; % Since we have Dirichlet boundary conditions, we don't need x=0
dX = 1/K; 

u = usolution(A,theta,lam,1,P0t(t)); 

lamt = 1+(1-u(end))/(2*(1-A(end)));  % compute lamt with the edges

tmpU = [ u; uf ] -lamt.*X ; % size K+1 x 1 


% add ghost node to A

tmpA = [2*D-A(1); A; 1]; % extend A by 1

FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;

F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;
F  = ( F(1:end-1) - F(2:end) ) / (dX); % F is the rhs to dA/dt

Alamrhs = F; 


Alamt = Alamrhs; 

% - solve

yt(1:K)= Alamt; 
% - solve

yt(K+1) = lamt ; 

yt = yt'; 

end

