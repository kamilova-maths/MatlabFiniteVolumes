% DATE:     2020 
% DESCR:    LargeOmega1AsymptoticsIsothermal
%           Main code looking at the fluctuations in the isothermal case (Gamma = 0).
%           Uses TimeDependentMOLIsothermal, ParametersDefinition, 
%           coupledPdeIsothermal, and usolution to solve pdes. 
%           First we calculate an analytical solution, where we choose the
%           Stokes number such that lambda  = 1. Then we calculate u for
%           this Stokes number and compare and print the difference. The
%           discrepancy is attributed to discretisation errors, and there
%           is no way to circumvent it. Then we use either simple or steady
%           conditions to start the code. We impose P0t as a function, but
%           can be chosen to be a constant value if required. The results
%           are plotted at discrete timesteps to show convergence towards
%           steady state. 
% INPUT: 
%           No input variables
%          
% OUTPUT:   Main outcomes: 
%           Acel: The N x 2*K matrix storing the cell values for A 
%           uinterf: The N x (2*K +1) matrix storing all interface values for
%           u
%           
% ADDITIONAL COMMENTS: 
%
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           coupledPdeIsothermal: This is where the isothermal pdes are, that will be then
%           solved with method of lines with ode15s


close all 



redo = 1;
while redo == 1
    TimeDependentMOLIsothermal % This gives Acel, temp, and uinterf, which have N rows and K columns
    prompt = ' Would you like to refine this further? (yes == 1) \n [Remember to rename SSNEW] ';
    redo = input(prompt);
end

global N K L P0 P0t    

% Our steady state should match the analytical steady state obtained in 
% TimeDependentMOLIsothermal , and we have to transpose them so that they
% are row vectors
A0 = A0';
u0bar  = 1./A0; 
     

% We define the sinusoidal Pin, with respect to the constant P0 used
% for the steady state
%% FOR SINUSOIDAL P


%tau      = linspace(0,2*pi,N)'; 
tau      = linspace(0,2*pi,N)';
Pin      = P0t(tau/omega);
Pinbar   = mean(P0t(tau/omega));
Pintilde = Pin- Pinbar;
Pintau   = Pintilde.*tau;

PinAn      = @(tau) P0+DeltaP*sin(tau); 
PintauAn   = @(tau) DeltaP*sin(tau).*tau; 
PintildeAn = @(tau) DeltaP*sin(tau);
PinbarAn   = @(tau) P0; 
% top , from 0 to lam0

xvec = linspace(0,1,K);
xtop = lam0.*xvec; 
dx   = (xtop(end)-xtop(1))/(K-1);
%bottom, from lam0 to 1
%xbot = lam0 + (L-lam0).*xvec; % for the integrals, this is okay, as lam0 is a constant.
% Integral from lam0 to x, where x goes from lam0 to 1 .The negative of
% this is the integral we want, i.e. integral from x to lam0 where x goes
% from lam0 to 1


phitau2  = zeros(N,1);
phitauAn = phitau2;

%tau = linspace(0,2*pi,N)';
dt = 1/(N-1);
% analytic expression for Pintau
term1An = (1/(2*pi))*integral(PintauAn,0,2*pi);

%vector expression for Pintau
term1   = (1/(2*pi))*trapz(Pintau).*((2*pi-0)*dt); 

for i = 1:N
    phitauAn(i) = term1An + integral(PintildeAn,0,tau(i)); 
    phitau2(i)  = term1 + trapz(Pintilde(1:i))*(tau(i)-0)*(1/(i)); 
end
phitau = term1 + cumtrapz(Pintilde)*(2*pi-0).*dt; 
check  = cumtrapz(Pintilde)*(2*pi-0).*dt; 

fA0     = flip(A0(1:K)); 
% The upper limit in this integral is 1 because dx already contains the
% lam0 factor, remember that xtop goes from 0 to lam0, so dx is lam0-0/K,
% and so the upper limit is 1*xtop(end)
Cint    = flip(-cumtrapz(1./(fA0)).*(-1)*dx);  
fac1    = smooth(derivative(A0(1:K).*Cint,dx));
fac2    = smooth(derivative(A0(1:K),dx).*Cint + A0(1:K).*derivative(Cint,dx)); 
A1tilde1 = -(1/3)*phitau*fac2';
A1tilde = (1/3)*phitau*(1-derivative(A0(1:K),dx).*Cint);

 
%A0prime = derivative(A0,dx); 
A0primelam0 = (St.*trapz([A0(1:K), 1]).*(1).*dx + Pinbar)./(3*(1)); 

%lam1tilde1 = -phitau*(1./(3*A0prime(K)*exp(-Gamma*th0(K))));
lam1tilde = -phitau*(1./(3*A0primelam0));


u0tilde = (Pintilde./3)*Cint(1:K);

tval = t(end)-2*pi/omega; 
indx = find(t>=tval,1);
indx = indx -1; 

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

xmatrix    = [lam*xcel',lam + (L-lam)*xcel'];
xmatrixint =  [lam*xint',lam + (L-lam)*xint(2:end)'];
xmat2      = lam*xcel';
xmat1int   = lam*xint';


Anum     = zeros(N,K);
u0barnum = zeros(N,K); 
for i = 1:N
    Anum(i,:) = interp1(xmatrix(i,:),Acel(i,:),linspace(0,L,K),'pchip'); 
end

xforu = linspace(0,1,K+1)*lam0;
for i = 1:N
     u0barnum(i,:) = interp1(xmatrixint(i,:),uinterf(i,:),linspace(0,L,K),'pchip'); 
end


omegat      = omega*t(indx:end);

AAvg        = mean(Anum(indx:end,:));

u0Avg       = mean(u0barnum(indx:end,:));
u0tildenum  = u0barnum - ones(N,1)*u0Avg; 


A1tildenum  = omega*(Anum-ones(N,1)*AAvg);

A1tildenum2 = smooth(A1tildenum);

lam1tildenum = omega*(lam-mean(lam(indx:end))); 


%% Calculating the analytical solution of the 'extra bit' for Acal


% x     = linspace(0,lam0,K);
F   = sqrt(6*St*D-P0^2);
a = atan(P0/F); 
Pstar =  6*a/(F); 
Aex   = @(x)(F^2/(6*St)).*tan(F.*(x+Pstar)./6).^2 - P0.^2./(6*St)+D ; 
uex = @(x) 1./Aex(x); 
c0 = -DeltaP; 



phiex = @(tau) -DeltaP*cos(tau) ; 
% Now we define the derivatives explicitly to prevent approximation errors
% due to finite differences

dAex = @(x) F^3.*(sec((F*x)./6 + a).^2).* tan((F.*x)/6 + a)./(18*St); 

ddAex = @(x)-((F^4 .*(-2 +cos((F.*x)./3 + 2.*a)).*(sec((F.*x)/6 + a).^4))./(108*St)); 


psiex = @(x) (-P0 + F.*tan((F.*x)./6 + a ))./St ;

uhatex = @(x) (St.*(F.*(-x + lam0) - 3*sin((F.*x)./3 + 2.*a) +  ...
    3.*sin((F.*lam0)./3 + 2.*a)))./F^3 ;

duhatex = @(x)(St.*(-F - F.*cos((F.*x)./3 + 2 .*a)))./F^3 ;
c1 = (uex(0)/uhatex(0))*(dAex(0)*uhatex(0) - (1/3)); 
Xhat =  @(x) omega.*psiex(x); 

Acal = @(x,tau) -((c1.*(Aex(x).*(Xhat(x).*c0.*D.*uhatex(0)) + ... 
    (c0 - phiex(tau))*(Aex(x).^2).*uhatex(x) + ... 
    c0.*psiex(x).*(1 + uhatex(x).*dAex(x))))./Aex(x)) ;



% We run the plotting code
run('PlottingFiles/ContoursLargeOmega1')


% Extra plot for Acal
figure; 

arg = @(tau) -tau - D*uhatex(0)*phiex(tau);
v0ex0 = -(D*St/(3*F^2))*(6-P0*lam0-(1/F)*(3*P0*sin((1/3)*F*lam0 + 2*atan(P0/F)))) ;

Fx = @(tau) v0ex0*phiex(tau);

plot(arg(tau),Fx(tau)); 
contourf(taumat, xmat1, Acal(xtop,tau), 100,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
xlim([0,2*pi])
ylim([0 lam0])
colorbar

if sav==1
    axis off
    colorbar off
    print(gcf, '-dpng', '-r300', '-painters', 'PlottingFiles/SavedPlots/AcalGzerou1zero.png')
end

return

i = 0;
delta = 9; 
while(i<N2)
    plot(xmat(indx+i,:), residue(i+1,:))
    hold on
    i = i+delta;
end
% Extra stuff that I should tidy-up / delete soon

u0barfun = @(x) interp1(xtop,u0bar,x);  
Cintfun = @(x) interp1(xtop,Cint,x);
psiint  = cumtrapz(1./(u0bar)).*dx;  

psifun = @(x) interp1(xtop,psiint,x); 

% integral from 0 to tau of u0(x,tau) is u0bar + Cint(x) integral of
% Pintilde over 3 

%xi =1; 
% Code that plots the characteristics of the leading order problem 
xval = lam0;
Numel = 20;
xi = linspace(0,5,Numel); 
eta = linspace(0,5,Numel); 
Pintildeint = zeros(size(xi));
psiprime = derivative(psifun(x),dx);
psiprimefun = @(x) interp1(xtop,psiprime,x);
figure;
for j = 1:length(eta)
    for i = 1:length(xi)
        Pintildeint(i) = trapz(PintildeAn(eta+xi(i)))*(xi(i)-0)/(Numel);
    end    
    
    Xhat = (psiprimefun(xval)).*(Pintildeint*Cintfun(xval) + xi.*u0barfun(xval));
    tau_char = eta(j) + xi; 
    plot(tau_char,Xhat)
    hold on
end
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$\tau$','Interpreter','latex')
ylabel('$\hat{X}$','Interpreter','latex')
title('Characteristics at $x = \lambda^{(0)}$, for $\xi \in (0,5)$ and $\eta \in (0,5)$', 'Interpreter','latex')
print(gcf, '-dpng', '-r300', '-painters', 'PlottingFiles/SavedPlots/Charxlam0.png')


% Plotting G1(0,s)
s = -tau - (dA0dx0/(Pinbar))*(Cint(1)*phitauAn - term1An);
G1x0 = (phitauAn/3)*(dA0dx0*Cint(1)-1);
plot(s,G1x0)
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$s$','Interpreter','latex')
ylabel('$G_1(0,s)$','Interpreter','latex')
print(gcf, '-dpng', '-r300', '-painters', 'PlottingFiles/SavedPlots/G1x0.png')


% Plotting G1(0,s) for analytical form of phitau for sinusoidal Pin
Nred = 500;
taured = linspace(0,2*pi,Nred)';
alpha = - (dA0dx0/(Pinbar))*Cint(1); 
beta = -alpha*term1;
s = -taured + alpha*(-DeltaP*cos(taured)) +beta;
%G1atx0 = -DeltaP*cos(-alpha*DeltaP*cos(taured*ones(1,Nred)) - beta - ones(Nred,1)*s)*(-(Pinbar/3)*alpha-1); 
G1atx0 = -DeltaP*cos(alpha*(-DeltaP*cos(taured)) +beta - s)*(-(Pinbar*alpha/3) -1/3);

figure;
contourf(taured*ones(1,Nred), ones(Nred,1)*s, G1atx0, 100,'LineColor', 'none')
% ax = gca;
% ax.YDir = 'reverse';
%xlim([0,2*pi])
%ylim([0 lam0])
colorbar


%Try to plot G_1(x,s)

dA0dx0 = Pinbar/(3*u0bar(1));
alpha = -(dA0dx0/(Pinbar))*Cint(1); 
beta = -alpha*term1;

psi = cumtrapz(1./u0bar)*dx ; 
hatX = omega*psi;
uhat0  = (1/3)*Cint;

s = ones(N,1)*hatX - tau*ones(1,K) -  (ones(N,1)*(uhat0./u0bar)).*(phitau*ones(1,K) - term1) ;

%c1s = phitau*Pinbar*(Cint(1)/9 - 1/(3*dA0dx(1)));

c1s = -DeltaP*cos(-alpha*DeltaP*cos(tau*ones(1,K)) + beta - s)*(-(1/3)-alpha*Pinbar/3)*u0bar(1);
rhs = c1s./u0bar; 

A1MS = rhs + A1tilde;
%A1tilde + (1/omega)*rhs
% ANALYTICAL
figure; 
contourf(taumat, xmat1, rhs , 100,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
xlim([0,2*pi])
ylim([0 lam0])
% if an == 1
%     if resc == 0
%         a = min(min(A1tilde));
%         b = max(max(A1tilde));
%     else
%         a = min(min(A1tildenum(indx:end,:)));
%         b = max(max(A1tildenum(indx:end,:))); 
%     end
% else
%         a = min(min(A1tilde));
%         b = max(max(A1tilde));
% end
% caxis([a, b])
% disp(['min A is ', num2str(a) , '  ', 'max A is ', num2str(b)])
colorbar