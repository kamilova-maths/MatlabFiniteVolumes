function dydx = odefunNH(x,y,alpha,Q,x1,x2,eps,St,Ta,Bi,Pe,lam)

% ODE function solver without the Heaviside. Here A can grow to be larger
% than one, and we are NOT imposing A=1 at the end. This should match the
% time dependent problem.
lam=1; 

P = y(1,:);
A = min(y(2,:),1);
T = y(4,:);
J = y(3,:);

% viscosity
Q = Q*(x>x1).*(x<x2);    % heat source
%He = (A<1);                     % Heaviside function
% He = 1/2*(1+tanh((1-A)/p.eps)); % regularised Heaviside function
%He = max(1-exp(-(1-A)/eps),0);% regularised Heaviside function
dydx = [ St*A.*lam^2;
        %(A.*P./(3*exp(-alpha*T)));
        (A.*P./(1.6));
        Pe*((J./A) - A.*Q) - 2*Bi*A.^(1/2).*(Ta-T);
        J./A
        ];
end