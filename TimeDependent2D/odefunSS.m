function dydx = odefunSS(x,y)
global Gamma x1 x2 Q tha Bi Pe eps2 St
P = y(1,:);
A = min(y(2,:),1);
T = y(4,:);
J = y(3,:);
mu = exp(-Gamma*T);  
% viscosity
Qfun = Q*(x>x1).*(x<x2);    % heat source
                     % air temperature
%He = (A<1);                     % Heaviside function
% He = 1/2*(1+tanh((1-A)/p.eps)); % regularised Heaviside function
He = max(1-exp(-(1-A)/eps2),0);% regularised Heaviside function
dydx = [ St*A;
        (A.*P./(3*mu)).*He;
        Pe*((J./A) - A.*Qfun) - 2*Bi*A.^(1/2).*(tha-T);
        J./A
        ];
end