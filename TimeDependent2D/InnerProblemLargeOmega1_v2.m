clear all
ParametersDefinition
global P0t Pinbar Pintilde phitau

%% Calculating P0 solution 
close all 

% We import the steady state of the OSCILLATING problem. This is different
% to the steady state version of the time dependent problem with P0
% constant. Change filename to match what we want to import: 
% Options for filenames: Regular parameter values with DeltaP = 0.8,
% SSDataOmega100All1.csv : Bi = Pe = P0 = Gamma = 1, omega = 100
% SSDataOmega300All1.csv : Bi = Pe = P0 = Gamma = 1, omega = 100

data = csvread('TextFiles/PeBiGamma1Omega100Qexp.csv');

A0 = data(2*K+1:3*K); 
th0 =  data(4*K+1:5*K);
phi0 = data(5*K+1:6*K); 
lam0 = data(10*K+3); 
A0lamt = A0.*lam0;

%% START FROM HERE WHEN YOU HAVE ALREADY CALCULATED STEADY STATE 

% Initial conditions 
y0(1:K) = A0lamt;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = phi0;		
y0(3*K+1) = lam0; 

%P0t = @(t)P0*cos(2*pi*t/0.1234);
%P0 = 1; 
omega = 100;
DeltaP = 0.8/P0; 
n = 100;
T= 2*pi*n/omega;
%T = 100; 

%DeltaP = 0.2;
P0t = @(t) P0 + DeltaP.*sin(omega*t); % base case 

%P0t = @(t) P0 + P0*sin(pi*t); % fewer oscillations than base case
%P0t = @(t) P0 + P0*sin(4*pi*t); % more oscillations than base case

% Independent variable for ODE integration 
tspan = [0 T];
tout = linspace(0,T,N);

%% ODE integration 
options = odeset('RelTol',1.0e-3,'AbsTol',1.0e-6);

tic
%[t,y] = ode15s(@coupledPde,tspan,y0); 
[t,y] = ode15s(@coupledPde,tout,y0,options); 
toc
N = length(t);
Alamt  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)
lam = y(:,3*K+1); 
A = Alamt./lam;
th = y(:,K+1:2*K)./A;
phi   = y(:,2*K+1:3*K);

% We calculate u with the solution for A, th and lam
u = zeros(size(A));
for i=1:N
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i),1,P0t(t(i)));   
end

% We add the Dirichlet boundary conditions 
thtop = [ zeros(N,1), ...
       th         ];

Acel = [ A, ones(N,K)];		% A (cell values)
Aint = ([D*ones(N,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)

uint  = [ u , ...
					uf.*ones(N,K+1) ];
temp = [th, phi];		% complete temperature profile (theta and phi)

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

plot(t,lam)
%return
st=1;    
if st==0
    xvector1 = [xcel*lam(end);lam(end) + xcel*(L-lam(end))];
    xvector2 = [xint*lam(end);lam(end) + xint(2:end)*(L-lam(end))];
    SS = [xvector1; Acel(end,:)'; temp(end,:)'; xvector2; uint(end,:)'; lam(end-1); P0t(t(end))]; 
    csvwrite('PeBiGamma1Omega100Qexp.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end


%% Calculating the Analytical oscillations, as well as the Inner Layer
% These are column vectors

    % Change filename to match what we want to import: 
    %Options for filenames: Regular parameter values with
    
    % SSDataP01.csv :  P0 = 1  (K=300)
    % SSK300P0p5.csv :  P0 = 0.5 (K=300)
    % SSKDataPeBi1.csv : P0 = 1, Pe = Bi = 1 (K=300)
    % SSDataP0Bi10.csv : P0 = 1, Bi = 10 (K=300)
    % SSDataP0Bi27.csv : P0 = 1, Pe = Bi = 27 (K=300)
    % SSDataP0Bi100.csv : P0 = 1, Bi = 100 (K=300)
    % SSDataP0BiPe10.csv : P0 = 1, Bi = Pe = 10 (K = 300) 
    % SSDataP01Gamma15.csv : P0 = 1, Gamma = 15 
    % SSDataPeBiGamma1.csv : P0 = Pe = Bi = Gamma = 1
    % SSDataPeBiGamma10.csv: P0 = Pe = Bi = 1, Gamma = 10
    % SSDataPeBiGamma1K500.csv : P0 = Pe = Bi = Gamma = 1, K = 500
    % SSDataPeBi1Gamma1K300Qexp.csv : P0 = Pe = Bi = Gamma = 1, K =300, Q
    % is a Gaussian (smooth continuous function of x)
    % SSDataPeBiB27Gamma23K300Qexp.csv : P0 = 1, Q is a Gaussian 
    
% Change filename to match what we want to import 
data = csvread('TextFiles/SSDataPeBi1Gamma1K300Qexp.csv');


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

PinAn =@(tau) P0+DeltaP*sin(tau); 
PintauAn = @(tau) DeltaP*sin(tau).*tau; 
PintildeAn = @(tau) DeltaP*sin(tau);
PinbarAn = @(tau) P0; 
% top , from 0 to lam0
dx = 1/(K-1);
xvec = linspace(0,1,K);
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
fac1 = derivative(A0(1:K).*Cint,dx);
fac2 = derivative(A0(1:K),dx).*Cint + A0(1:K).*derivative(Cint,dx); 
A1tilde = -(1/3)*phitau*fac1;

dth0dx = smooth(derivative(th0(1:K),dx))';
t1tilde = -(1/3).*phitau*Cint.*(dth0dx);

%A0prime = derivative(A0,dx); 
A0primelam0 = (St.*trapz([A0(1:K), 1]).*(lam0).*dx + Pinbar)./(3*exp(-Gamma*th0(K))); 

%lam1tilde1 = -phitau*(1./(3*A0prime(K)*exp(-Gamma*th0(K))));
lam1tilde = -phitau*(1./(3*A0primelam0*exp(-Gamma*th0(K))));


u0tilde = (Pintilde./3)*Cint(1:K);

tval = t(end)-2*pi/omega; 
indx = find(t>=tval,1);
indx = indx -1; 

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

xmatrix = [lam*xcel',lam + (L-lam)*xcel'];
xmatrixint =  [lam*xint',lam + (L-lam)*xint(2:end)'];

xmat1int = lam*xint';


Anum = zeros(N,K);
thnum = zeros(N,K);
u0barnum = zeros(N,K);
%xtop = lam0.*linspace(xint(2)/2,1-xint(2)/2,K); 
xtop = lam0.*linspace(0,1,K); 
for i = 1:N
    Anum(i,:) = interp1(xmatrix(i,:),Acel(i,:),linspace(0,L,K),'spline'); 
    thnum(i,:) = interp1(xmatrix(i,:),temp(i,:),linspace(0,L,K),'spline');
end

xforu = linspace(0,1,K+1)*lam0;
for i = 1:N
     u0barnum(i,:) = interp1(xmatrixint(i,:),uint(i,:),linspace(0,L,K),'spline'); 
end
 
omegat = omega*t(indx:end);
AAvg = mean(Anum(indx:end,:));
thAvg = mean(thnum(indx:end,:));
u0Avg = mean(u0barnum(indx:end,:));
u0tildenum = u0barnum - ones(N,1)*u0Avg;

A1tildenum = omega*(Anum-ones(N,1)*AAvg);

t1tildenum = omega*(thnum-ones(N,1)*thAvg);  

lam1tildenum = omega*(lam-mean(lam(indx:end))); 

sav = 0;
close all
%plot(t,lam)
%PlottingLargeOmega1Asymptotics

X       = xtop.*sqrt(omega); 
dA0dx   = derivative(A0,dx);

ddA0dx  = derivative(dA0dx,dx); 
%dth0dx  = derivative(th0,dx); 
ddth0dx = derivative(dth0dx,dx);
cA0     = D.*ones(size(X));
cA1     = dA0dx(1).*X; 
cA2     = (1/2).*(X.^2).*ddA0dx(1) +(1/3).*phitau*(1-dA0dx(1).*Cint(1));
u0hat0  = (1/3).*Cint(1);

u0xxhat0 = (1/(3*D)).*(dA0dx(1)/D - Gamma.*dth0dx(1));
u0xxbar0 = (Pinbar/(3*D)).*(dA0dx(1)/D - Gamma.*dth0dx(1)) - St/3; 
U0 = Pinbar./(3*dA0dx(1)) + u0hat0.*Pintilde;
U1 = -(Pin./(3*D))*X; 
U2x = X.*u0xxbar0 + Pintilde*(u0xxhat0.*X);
th1hat0 = -(dth0dx(1)/3).*Cint(1); 

cA3dt = -U1.*dA0dx(1) - (-Pin./(3*D))*cA1 - U0*(X.*ddA0dx(1)) - U2x.*(ones(N,1)*cA0);

cA3 = zeros(N,K); 
for i = 1:K
    cA3(:,i) = cumtrapz(cA3dt(:,i))*(2*pi-0)*dt;
end

opt = 'red2';
switch opt
    case 'red1'
        %% Reduction 1

        %t1initialT = (DeltaP*th1hat0)/(-X(end)).*(X) + th1hat0*DeltaP; 
        t1initialTalt = th1hat0*DeltaP*exp(-X); 
        %y0T = t1initialT(2:end-1); % I don't solve at the edge points as we have Dirichlet there 
        y0T = t1initialTalt(2:end-1); % I don't solve at the edge points as we have Dirichlet there 
        k = 5; % to get it to converge to a final solution. 

        tout2 = linspace(0,k*2*pi,N);
        dA0dx0= D*Pinbar/3;
        Qfun = @(x) Q*(x>x1).*(x<x2);    % heat source  
        Qfunvec = Qfun(xtop);
        %dth0dx0 =  (2*Bi/D).*trapz(sqrt(A0).*(th0(1:K)-tha) - Qfunvec.*A0).*(lam0).*dx;
        
        dth0dx0 = ((1/Pe).*D*ddth0dx(1) - (2*Bi/Pe).*sqrt(D)*(-tha) )./(1-(1/Pe).*dA0dx0);
        
        [t2,yTheta] = ode15s(@(t,y)InnerLayerT(t,y,dth0dx0,dA0dx0,u0hat0,th1hat0,ddth0dx),tout2,y0T); 

        N2 = length(t2); 
        phitauN2 = DeltaP.*cos(t2); 
        %phitauN2 = interp1(tau,phitau,t2); % This includes the DeltaP factor

        yTheta = [phitauN2.*th1hat0, yTheta, zeros(N2,1) ];
        indxa = find(t2>=(k-1)*2*pi,1);
        indxa = indxa -1 ; 

        % Readjusting the sizes of things so that we can plot the composite
         N2 = length(t2(indxa:end)); 

        taucomp = linspace(0,2*pi,N2)'; 

        Pin = P0t(taucomp/omega);
        Pinbar2 = mean(P0t(taucomp/omega));
        Pintilde2 = Pin- Pinbar2;
        Pintau2 = Pintilde2.*taucomp;
        dt2 = 1/N2;
        term12 = (1/(2*pi))*trapz(Pintau2).*((2*pi-0)*dt2); 
        phitau2= term12 + cumtrapz(Pintilde2)*(2*pi-0).*dt2; 

        t1tilde2 = -(1/3).*phitau2*Cint.*(dth0dx);
        composite1 = t1tilde2 + yTheta(indxa:end,:); 

        xmat2 = ones(N2,1)*linspace(0,lam0,K); 
        Xmat = ones(N2,1)*X; 
        tau2 = linspace(0,2*pi,N2)';
        taumat2 = tau2*ones(1,K); 

        % % PLOTTING COMPOSITE
         figure; 

       %contourf(taumat2, xmat2, yTheta(indxa:end,:), 100,'LineColor', 'none')
        contourf(taumat2, xmat2, composite3, 100,'LineColor', 'none')

        ax = gca;
        ax.YDir = 'reverse';
        xlim([0,2*pi])
        ylim([0 lam0])
        %title('$\tilde{A}_1$ Analytic, $\omega_1=10$','Interpreter','latex','FontSize',14);
        %caxis([min(min(A1tilde)),max(max(A1tilde))])
        %caxis([min(min(A1tildenum)),max(max(A1tildenum))])
        colorbar
        %caxis([min(min(yTheta(indxb:end,:))),max(max(yTheta(indxb:end,:)))])
        %caxis([min(min(t1tildenum)),max(max(t1tildenum))])



    case 'red2'
        %% Reduction 2


        t1initialThat = zeros(size(X));
        y0T = t1initialThat(2:end-1); % I don't solve at the edge points as we have Dirichlet there 

        % The solution to this problem is That, and Ttilde = That + W, W given
        % below. Then That is the solution to purple PDE, where Pin is sinusoidal
        % (so this only works for this case) 
        k = 5; % to get it to converge to a final solution. 
        tspan2 = [0 k*2*pi];
        tout2 = linspace(0,k*2*pi,N);
        dA0dx0= D*Pinbar/3;
        Qfun = @(x) Q*(x>x1).*(x<x2);    % heat source  
        Qfunvec = Qfun(xtop);
        %dth0dx0 =  (2*Bi/D).*trapz(sqrt(A0).*(th0(1:K)-tha) - Qfunvec.*A0).*(lam0).*dx;
        
        dth0dx0 = ((1/Pe).*D*ddth0dx(1) - (2*Bi/Pe).*sqrt(D)*(-tha) )./(1-(1/Pe).*dA0dx0);
        
        %[t2,yThetaThat] = ode15s(@(t,y)InnerLayerThat(t,y,X(2:end-1)'),tspan2,y0T);
        [t2b,yThetaThat] = ode15s(@(t,y)InnerLayerThat(t,y,X(2:end-1)',dth0dx0,dA0dx0,u0hat0,th1hat0,ddth0dx,DeltaP),tout2,y0T); 

        N2b1 = length(t2b); 

        yThetaThat = [zeros(N2b1,1), yThetaThat, zeros(N2b1,1)];

        W = DeltaP.*(th1hat0.*cos(t2b))*exp(-X);

        % Now yTheta is 
        yThetaT = yThetaThat + W;
        indxb = find(t2b>=(k-1)*2*pi,1);
        indxb = indxb -1 ; 

        % Readjusting the sizes of things so that we can plot the composite
         N2b2 = length(t2b(indxb:end)); 
        taucomp = linspace(0,2*pi,N2b2)'; 
        Pin = P0t(taucomp/omega);
        Pinbar2 = mean(P0t(taucomp/omega));
        Pintilde2 = Pin- Pinbar2;
        Pintau2 = Pintilde2.*taucomp;
        dt2 = 1/N2b2;
        term12 = (1/(2*pi))*trapz(Pintau2).*((2*pi-0)*dt2); 
        phitau2= term12 + cumtrapz(Pintilde2)*(2*pi-0).*dt2; 

        t1tilde2 = -(1/3).*phitau2*Cint.*(dth0dx);
        % 
        composite2 = t1tilde2 + yThetaT(indxb:end,:); 

        %composite2resc = abs(min(min(t1tildenum))).*composite2./abs((min(min(composite2))));
        % composite2 = t1tilde2 + yTheta(indx2:end,:); 

        xmat2 = ones(N2b2,1)*xtop; 
        %Xmat = ones(N2,1)*X; 
        tau2 = linspace(0,2*pi,N2b2)';
        taumat2 = tau2*ones(1,K); 
        T1 = @(X,tau) real((DeltaP.*th1hat0(1).*cos(tau).*exp(-sqrt(1i*Pe).*X)));

        composite3 = t1tilde2 + T1(xmat2.*sqrt(omega),taumat2);
        % % PLOTTING COMPOSITE
         figure; 

        %contourf(taumat2, xmat2, yThetaT(indxb:end,:), 100,'LineColor', 'none')
        contourf(taumat2, xmat2, composite3, 100,'LineColor', 'none')
        ax = gca;
        ax.YDir = 'reverse';
        xlim([0,2*pi])
        ylim([0 lam0])
        colorbar
        caxis([-0.56912, 0.54405])
        if sav==1
            axis off
            colorbar off
            print(gcf, '-dpng', '-r300', '-painters', 'PlottingFiles/SavedPlots/t1tildeComposite.png')
   
        end
        %caxis([min(min(t1tildenum)),max(max(t1tildenum))])
        %caxis([min(min(yTheta(indxb:end,:))),max(max(yTheta(indxb:end,:)))])

end

%T1 = @(X,t) -(DeltaP.*th1hat0(1).*cos(t).*exp(-sqrt(Pe).*X));
        % % PLOTTING COMPOSITE
         figure; 

        %contourf(taumat2, xmat2, T1(xmat2.*sqrt(omega),taumat2), 100,'LineColor', 'none')
        contourf(taumat2, xmat2, composite3, 100,'LineColor', 'none')
        ax = gca;
        ax.YDir = 'reverse';
        xlim([0,2*pi])
        ylim([0 lam0])
        colorbar
        if sav==1
            axis off
            colorbar off
            print(gcf, '-dpng', '-r300', '-painters', 'SavedPlots/t1tildeAnalytic.png')
   
        end
        
        % for a particular form of Pintau, we have the following solution
        % for A 
 U0hat = @(tau) Pinbar./(3*dA0dx0) + (1/3).*Cint(1).*DeltaP.*sin(tau);        
A2inner = @(X,tau) exp(-X).*((-Pinbar - 1 )./(3*U0hat(tau))).*X - (1/3).*(1-dA0dx0.*sin((tau-1./U0hat(tau)).*X)).*Cint(1); 



%composite3 = t1tilde + T1((1/sqrt(omega)).*X,tau);
