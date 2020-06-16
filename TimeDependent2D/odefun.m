function dydx = odefun(x,y,eps)

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