function dydx = odefunNH(x,y)

% DATE:  2020

% DESCR: Function that contains the ODEs for the steady state problem for 
%        coupled flow and temperature problem without Heaviside. 
%        We use global variables. It is in the required format for being 
%        used with bvp4c and solinit
%
% INPUT: x: vector of x values 
%        y: matrix contaiting P, A, T, and J (in that order)
%           
% OUTPUT: dydx: Matrix with dydx values at specific x.    

% ADDITIONAL COMMENTS: 
%           

% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           bvpinit, bvp4c: to solve the ODEs
%           InitialConditionsSteady: Main program for calculating steady
%           state solutions.
%           odefun: Contains the equations to solve, for H == 1
%           bcfun: Contains the boundary conditions for the ODEs. 

global Pe Bi tha Gamma St x1 x2 Q 
P = y(1,:);
A = min(y(2,:),1);
T = y(4,:);
J = y(3,:);

mu = 3*exp(-Gamma*T);

% viscosity
Q = Q*(x>x1).*(x<x2);    % heat source
dydx = [ St*A;
        (A.*P./(mu));
        Pe*((J./A) - A.*Q) - 2*Bi*A.^(1/2).*(tha-T);
        J./A
        ];
end