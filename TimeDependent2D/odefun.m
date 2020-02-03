function dydx = odefun(x,y,alpha,Q,x1,x2,eps,St,Ta,Bi,Pe)
P = y(1,:);
A = min(y(2,:),1);
T = y(4,:);
J = y(3,:);
mu = exp(-alpha*T);  
% viscosity
Q = Q*(x>x1).*(x<x2);    % heat source
%He = (A<1);                     % Heaviside function
% He = 1/2*(1+tanh((1-A)/p.eps)); % regularised Heaviside function
He = max(1-exp(-(1-A)/eps),0);% regularised Heaviside function
dydx = [ St*A;
        (A.*P./(3*mu)).*He;
        Pe*((J./A) - A.*Q) - 2*Bi*A.^(1/2).*(Ta-T);
        J./A
        ];
end