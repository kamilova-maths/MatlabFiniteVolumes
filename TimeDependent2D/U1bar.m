% Code to determine u1bar. We have one first order and one second order ODE
% to solve for u1bar and A1bar, for the case where A1bar includes the
% multiple scales part. This is the isothermal case, and we need to already
% have the solutions u0bar, A0, dA0dx, du0dx, Pinbar, St, and all of that lovely
% stuff%

% We want to find what condition for U1bar at x =0 is required such that
% U1bar(lambda0) is zero .

% We first define the analytical expressions for the solution of the
% isothermal problem 

ParametersDefinition
DeltaP = 0.8;
%St = 1; 
lamex = @(St) atan(sqrt(P0^2 + 6*St - 6*D*St)/sqrt(6*St*D-P0^2)).*6./(sqrt(6*St*D-P0^2)) - ...
        6*atan(P0./sqrt(6*St*D-P0^2))./(sqrt(6*St*D-P0^2)); 

lam0 = lamex(St);% We use St = 1, even though this means that lam0 > 1

% x     = linspace(0,lam0,K);
F   = sqrt(6*St*D-P0^2);
a = atan(P0/F); 
Pstar =  6*a/(F); 
Aex   = @(x)(F^2/(6*St)).*tan(F.*(x+Pstar)./6).^2 - P0.^2./(6*St)+D ; 
uex = @(x) 1./Aex(x); 

% Now we define the derivatives explicitly to prevent approximation errors
% due to finite differences

dAex = @(x) F^3.*(sec((F*x)./6 + a).^2).* tan((F.*x)/6 + a)./(18*St); 

ddAex = @(x)-((F^4 .*(-2 +cos((F.*x)./3 + 2.*a)).*(sec((F.*x)/6 + a).^4))./(108*St)); 


psiex = @(x) (-P0 + F.*tan((F.*x)./6 + a ))./St ;

uhatex = @(x) (St.*(F.*(-x + lam0) - 3*sin((F.*x)./3 + 2.*a) +  ...
    3.*sin((F.*lam0)./3 + 2.*a)))./F^3 ;

duhatex = @(x)(St.*(-F - F.*cos((F.*x)./3 + 2 .*a)))./F^3 ;


% We define constants that we will need for the next steps 
%c0 = -1 ; % I am so sorry about this ... 
c0 = -DeltaP; 

c1 = (uex(0)/uhatex(0))*(dAex(0)*uhatex(0) - (1/3)); 
c1part = (dAex(0)*uhatex(0) - (1/3));
% We continue with the functions necessary for the coefficients in the
% U1bar problem 

% g1ex =@(x) (c1.*(Aex(x).^3).*psiex(x))./(Aex(x).^2 - psiex(x).*dAex(x)); 
% 
% dg1ex = @(x) c1.*Aex(x).^2 .*((Aex(x).^2 ).* psiex(x) .* dAex(x) -  ...
%     3 .* (psiex(x).^2).* (dAex(x).^2 )+ (Aex(x).^3).* Aex(x) +Aex(x) .* ...
%     (psiex(x).^2).* ddAex(x))./(Aex(x).^2 - psiex(x).* dAex(x)).^2 ;
% 
% 
% g2ex = @(x)  -((c1.*Aex(x).*(Aex(x).* (pi - c0.* D .*uhatex(0)) + ...
%     c0.*(Aex(x).^2).*uhatex(x) + c0 .*psiex(x).* ...
%     (1 + uhatex(x).* dAex(x))))./(Aex(x).^2 - psiex(x).*dAex(x)));
% 
% dg2ex = @(x) (c1./(Aex(x).^2 - psiex(x).* dAex(x)).^2).*(-c0.*(Aex(x).^4).*uhatex(x).*dAex(x) + ...
%     c0 .* (psiex(x).^2) .* (dAex(x).^2).* (1 + uhatex(x) .* dAex(x)) - ...
%     c0 .* (Aex(x).^5).* duhatex(x) + ...
%     Aex(x) .* psiex(x) .* ((dAex(x).^2).* (2 *pi - 2.*c0.*D.*uhatex(0) +...
%     c0 .* psiex(x) .* duhatex(x)) - c0 .* psiex(x) .* ddAex(x))+ ... 
%     (Aex(x).^2).* (dAex(x).*(psiex(x).* (c0 + 4.* c0 .*uhatex(x) .* dAex(x)) + (-pi + ...
%     c0 .* D .* uhatex(0)) .* Aex(x)) + (-pi + c0 .* D .* uhatex(0)) .* psiex(x) .* ddAex(x)) - ...
%     c0 .* (Aex(x).^3).* ((1 + 2 .*uhatex(x).* dAex(x)).* Aex(x) +2.*uhatex(x) .* psiex(x) .* ddAex(x))); 

%dA0dx = derivative(A0,dx);

% We define the coefficients for the ODE in u1, as well as partially define
% the RHS

% Alternative, if G1(x) = 0, then g1ex, dg1ex, g2ex, dg2ex are zero. 

g1ex =@(x) 0;
g2ex = @(x) 0;
dg1ex =@(x) 0;
dg2ex = @(x) 0; 

alpha1 = @(x) c1 + (2.*dAex(x))./Aex(x) + (c1.*(Aex(x).^2))./(-(Aex(x).^2) + psiex(x).*dAex(x)); 

alpha2 = @(x) (1/3 ).*((3.*dAex(x).*(2.*Aex(x).*dAex(x) - dg1ex(x)))./ ...
   Aex(x).^3 - (St.*(Aex(x).^3 )+ 6.* (dAex(x)).^2 - 3 .*Aex(x).* ddAex(x) )./(Aex(x).^2) + ( ...
   c1.*psiex(x) .* (St.*(Aex(x)).^3 + 6.*(dAex(x).^2) -3.*Aex(x).*ddAex(x)))./( (Aex(x).^3) - Aex(x) .*psiex(x).* dAex(x)));

beta1 = @(x) -((St.*g2ex(x))./(3.* Aex(x))) - (2.*g2ex(x).*dAex(x).^2)./(Aex(x).^4)+ (dAex(x) .*dg2ex(x))./(Aex(x).^3) + (g2ex(x).*ddAex(x)./(Aex(x).^3)); 
 
%beta2 = @(x) -((D.*St)./3) -(D.* dAex(x).^2)./(Aex(x).^3) +
%(D*ddAex(x)./(Aex(x).^2)); % beta2 is zero!
%u1atx0 = (D*pi + c0 *psiex(0) + c0*uhatex(0)*psiex(0)*dAex(0))/(D^2 *psiex(0));

v0 = (P0/(3*(D^2)))*g2ex(0); 

u1bar0 = fsolve(@(u1bar0)U1barODE(u1bar0, v0, alpha1, alpha2, beta1, D, uhatex(0), c0, c1, psiex(lam0), dAex(lam0), lam0), 0.01);

[xout, y] =ode45(@(x,y)odefunU1bar(x,y,alpha1, alpha2, beta1),[0 lam0],[u1bar0; v0]);

u1 = y(:,1);
v = y(:,2); 

figure
plot(xout, u1)


phiex = @(tau) -DeltaP*cos(tau) ; 

% We define G1(x) as a vector 

G1 = -((c1.*Aex(xout).*(Aex(xout) .*(pi - c0.*D.*uhatex(0)) + (Aex(xout).^2).*(c0.*uhatex(xout) - u1.*psiex(xout)) + ...
    c0.*psiex(xout).*(1 + uhatex(xout).*dAex(xout))))./((Aex(xout).^2) - psiex(xout).* dAex(xout))) ;
G1 = 0; 

% Plotting curly A = h1(x) -c1 Xhat + c1 uhat/u0 phitau
u1 =0; 
f = (c1.*Aex(xout).*(Aex(xout).*(c0.*D.* uhatex(0)) - c0 .*psiex(xout) + (Aex(xout).^2).* (uhatex(xout) .*(-c0) + u1 .*psiex(xout))) + (G1 - ...
    c0.*c1.*Aex(xout).*uhatex(xout)).*psiex(xout).*dAex(xout))./(Aex(xout).^2) ;
% Matrix, length(xout) rows and length(tau) columns
tau = linspace(0,2*pi,100);% column vector
omega = 100;
Xhat =  @(x) omega.*psiex(x); 
%Acal = f - c1.*Xhat + c1.*phiex(tau)*(uhatex(xout)./uex(xout));


Acal = @(x,tau) -((c1.*(Aex(x).*(Xhat(x).*c0.*D.*uhatex(0)) + ... 
    (Aex(x).^2).*uhatex(x)*(c0 - phiex(tau)) + ... 
    c0.*psiex(x).*(1 + uhatex(x).*dAex(x))))./Aex(x)) ;

contourf(tau,xout,Acal',50,'LineColor','none')
ax = gca;
ax.YDir = 'reverse';
xlim([0,2*pi])
ylim([0 lam0-1])
colorbar


function X = U1barODE(u1bar0, v0, alpha1, alpha2, beta1, D, uhatex0, c0, c1, psilam, dAexlam, lam0)  
%dx = (xtop(end)-xtop(1))./(length(xtop)-1);
%dA0dx = derivative(A0,dx);
%du0dx = derivative(u0bar,dx);
%du0dxdx = derivative(du0dx,dx); % there should be better ways to calculate these  
%du0dxdx = (-St*A0 - 3*dA0dx.*du0dx)./(3*A0);

% We write down the coefficients for the ODE


[ ~,y] = ode45(@(x,y)odefunU1bar(x,y,alpha1, alpha2, beta1),[0 lam0],[u1bar0; v0]);
%P1bar = y(:,2); 
u1bar = y(:,1);
u1bar0 = u1bar(1); 
% cond = D*(c0*c1*uhatex0 + u1bar0) - psilam *(c0*c1 +(c1 *pi + ...
%           D*u1bar0) *  dAexlam)/(1 - c1 + ...
%           dAexlam + psilam*(-c1 + (-1 + c1 - ...
%           dAexlam)*dAexlam)) ;
% there could be a minus sign missing here       
cond =  -(dAexlam.*(-D.*(c0.*c1.*uhatex0 + ... 
      u1bar0) + psilam.*(c0.*c1 + (c1.* pi + ...
         D.*u1bar0).*dAexlam)))./(c1 + ... 
 2*dAexlam.*(-1 + psilam.*dAexlam)) ; 

% cond = (dAexlam.*(-D.*(c0.*c1.*uhatex0 + ... 
%       u1bar0) + psilam.*((1 + c0).* c1 + (c1.*pi + ... 
%          D.*u1bar0- c1.*psilam).*dAexlam)))./(c1 +... 
%  2.*dAexlam.*(-1 + psilam.*dAexlam )) ; 

%cond = (D/2)*u1bar0; 
u1lam0 = u1bar(end); 
X = abs(u1lam0 - cond);  
% figure;
% plot(xout,A1bar)
% hold on 
% plot(xout,u1bar)

end


