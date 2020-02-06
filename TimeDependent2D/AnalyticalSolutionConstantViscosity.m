% Function to use fsolve to compute numerically what the integration
% constants should be in order to satisfy the boundary conditions. These
% equations are the evaluations of the analytical solution for P and A for
% constant viscosity at x=0, and equated to their boundary conditions. 

% Here, F(1) =P and F(2) = A, whereas c(1)=c1 and c(2)=c2, as defined in
% Mathematica printout

function F = AnalyticalSolutionConstantViscosity(c,St,P0,D)

F(1) =  sqrt(2).* sqrt (c(1)*(tanh(0.5.*sqrt((1/3).*c(1)).*c(2))).^2)-P0; 
F(2) = St.*D+(1/3).*c(1).*(sech(0.5.*(sqrt((1/3).*c(1)).*c(2)))).^2;

end