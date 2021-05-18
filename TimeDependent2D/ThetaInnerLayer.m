P0experiments % This gives Acel, temp, and uint, which have N rows and K columns
global N K Gamma L P0 P0t    
% These are column vectors

% Change filename to matth what we want to import 
%data = csvread('SSDataP01Bi27.csv');
A0 = data(2*K+1:3*K)'; 
u0bar = data(8*K+2:10*K+2)';
th0 =  data(4*K+1:6*K)';
    
lam0 = data(10*K+3); 

% We define the sinusoidal Pin, with respect to the constant P0 used
% for the steady state
%% FOR SINUSOIDAL P
tau = linspace(0,2*pi,N)'; 
Pin = P0t(tau/omega);
Pinbar = mean(P0t(tau/omega));
Pintilde = Pin- Pinbar;
Pintau = Pintilde.*tau;
%DeltaP = 0.5;  % not sure about this
%omegaP = 5;
PinAn =@(tau) P0+DeltaP*sin(tau); 
PintauAn = @(tau) DeltaP*sin(tau).*tau; 
PintildeAn = @(tau) DeltaP*sin(tau);
PinbarAn = @(tau) P0; 
% top , from 0 to lam0
dx = 1/(K-1);
xvec = linspace(0,1,K);
xtop = lam0.*xvec; 
%bottom, from lam0 to 1
xbot = lam0 + (L-lam0).*xvec; % for the integrals, this is okay, as lam0 is a constant.
% Integral from lam0 to x, where x goes from lam0 to 1 .The negative of
% this is the integral we want, i.e. integral from x to lam0 where x goes
% from lam0 to 1

%fac2 = derivative(A0.*Cint,dx); 
%tau = omegaP.*t;

phitau2 = zeros(N,1);
phitauAn = phitau2;

%tau = linspace(0,2*pi,N)';
dt = 1/(N-1);
% analytic expression for Pintau
term1An = (1/(2*pi))*integral(PintauAn,0,2*pi);

%vector expression for Pintau
term1 = (1/(2*pi))*trapz(Pintau).*((2*pi-0)*dt); 

for i = 1:N
 phitauAn(i) = term1An + integral(PintildeAn,0,tau(i)); 
 phitau2(i) = term1 + trapz(Pintilde(1:i))*(tau(i)-0)*(1/(i)); 
end
phitau= term1 + cumtrapz(Pintilde)*(2*pi-0).*dt; 
fth0 = flip(th0(1:K));
fA0 = flip(A0(1:K)); 
Cint = flip(-cumtrapz(1./(exp(-Gamma*fth0).*fA0)).*(-lam0)*dx); 
%fac1 = derivative(A0(1:K).*Cint,dx);
%A1tilde = -(1/3)*phitau*fac1;

dth0dx = smooth(derivative(th0(1:K),dx))';
t1tilde = -(1/3).*phitau*Cint.*(dth0dx);

%A0prime = derivative(A0,dx); 
%A0primelam0 = (St.*trapz([A0(1:K), 1]).*(lam0).*dx + Pinbar)./(3*exp(-Gamma*th0(K))); 

%lam1tilde1 = -phitau*(1./(3*A0prime(K)*exp(-Gamma*th0(K))));
%lam1tilde = -phitau*(1./(3*A0primelam0*exp(-Gamma*th0(K))));


%u0tilde = (Pintilde./3)*Cint(1:K);
% OPTION 2: Use the timedependent problem to find the steady state
t1initial = interp1(linspace(0,lam0,K),t1tilde(1,:),linspace(0,lam0,K-2));
y0 = -t1initial;
tspan = [0 2*pi];
[t,y] = ode15s(@(t,y)InnerLayer(t,y,t1tilde(:,1)),tspan,y0); 


N2 = length(t); 
t1tildefin = interp1(linspace(0,2*pi,N), t1tilde(:,1),t);
T1tilde = [zeros(N2,1), y, t1tildefin]; 

% tmatrix = t*ones(1,K);
% xmatrix = ones(N2,1)*linspace(0,1,K);
% contourf(tmatrix, xmatrix,T1tilde,20,'LineColor', 'none')
% ax = gca;
% ax.YDir = 'reverse';
% caxis([0 1])
% xlim([0 pi/10])

% Resizing the vectors to plot the composite solution 
xtemp = linspace(0,lam0,K); 
tau2 = linspace(0,2*pi,N2);
t1tilderesized = interp2(xtemp,tau,t1tilde,xtemp,tau2');


figure;
thcomposite = T1tilde + t1tilderesized - [zeros(1,K); t1tildefin(2:end)*ones(1,K)]; 
contourf(tmatrix, xmatrix,thcomposite,20,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
%caxis([0 1])
xlim([0 2*pi])



