% Code to determine u1bar inside an fsolve. The issue is that we know that
% P(0) = 0, but we only have U1bar at the other end (x=lambda0)
% We have one first order and one second order ODE
% to solve for u1bar and A1bar, for the case where A1bar includes the
% multiple scales part. This is the isothermal case, and we need to already
% have the solutions u0bar, A0, dA0dx, du0dx, Pinbar, St, and all of that lovely
% stuff%

function dydx = odefunU1bar(x,y,alpha1, alpha2, beta1)

% the order is y(1) = u1bar, y(2) = v
u1bar = y(1);
v = y(2); 


dydx = [ v;
        -alpha1(x)*v-alpha2(x)*u1bar-beta1(x)];
end
