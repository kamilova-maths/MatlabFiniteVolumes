function dydx = odefun(x,y,eps)
% DATE:  2020

% DESCR: Function that contains the ODEs for the steady state problem for 
%        coupled flow and temperature problem. We use global variables and
%        the value of eps provided by InitialConditionsSteady. It is in the
%        required format for being used with bvp4c and solinit
%
% INPUT: x: vector of x values 
 %       y: matrix contaiting P, A, T, and J (in that order)
 %       eps: determines how the heaviside transitions from 0 to 1
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
%           odefunNH: Contains the equations to solve, for H \neq 1
%           bcfun: Contains the boundary conditions for the ODEs. 

global Pe Bi tha Gamma St x1 x2 Q 
P = y(1,:);
A = min(y(2,:),1);
T = y(4,:);
J = y(3,:);
mu = exp(-Gamma*T);  
% viscosity
Qfun = @(x) Q*(x>x1).*(x<x2);    % heat source
%He = (A<1);                     % Heaviside function
%He = 1/2*(1+tanh((1-A)/eps)); % regularised Heaviside function
He = max(1-exp(-(1-A)/eps),0);% regularised Heaviside function
dydx = [ St*A;
        (A.*P./(3*mu)).*He;
        Pe*((J./A) - A.*Qfun(x)) + 2*Bi*A.^(1/2).*(T-tha);
        J./A
        ];
end