function res = bcfun(ya,yb,Tin)
% DATE:  2020
%
% DESCR:    bcfun(ya,yb,Tin)
%           Simple code that includes the boundary conditions for the
%           steady state problem for coupled flow and temperature problem. 
%
% INPUT: 
%           ya: Conditions at x = 0 (initial part of the domain)
%           yb: Conditions at x = L (last point of the domain)
%           Tin: Initial condition for thet a
% OUTPUT:   res: residue matrix for the bvp4c routine. This is how the 
%           boundary conditions are imposed in the steady state problem for
%           coupled flow and temperature.
% ADDITIONAL COMMENTS: 
%           There is a twin version of this code in the steady
%           state folder. Leave both there. I think it's worth having a 
%           functional example in the Time Dependent part as well. 
%
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           bvpinit, bvp4c: to solve the ODEs
%           odefun: Contains the equations to solve, for H == 1
%           odefunNH: Contains the equations to solve, for H \neq 1
%           bcfun: Contains the boundary conditions for the ODEs. 



global P0 D
res = [ ya(1,:)-P0(1);
        ya(2,:)-D; 
        yb(3,:);
        ya(4,:)-Tin 
        ];
end