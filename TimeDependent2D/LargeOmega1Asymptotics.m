% DATE:     2020 
% DESCR:    LargeOmega1Asymptotics
%           Main code looking at the fluctuations for large omega1 .
%           Uses TimeDependentMOL, ParametersDefinition, 
%           coupledPde, and usolution to solve pdes. 
%           We use either simple or steady
%           conditions to start the code. We impose P0t as a function, but
%           can be chosen to be a constant value if required. We plot the
%           contours with an external routine
% INPUT: 
%           No input variables
%          
% OUTPUT:   Main outcomes: 
%           Acel: The N x 2*K matrix storing the cell values for A 
%           uinterf: The N x (2*K +1) matrix storing all interface values for
%           u
%           theta: N x 2*K matrix storing cell values for theta
%           
% ADDITIONAL COMMENTS: 
%
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           coupledPde: This is where the des are, that will be then
%           solved with method of lines with ode15s



close all 

redo = 1;
while redo == 1
    TimeDependentMOL % This gives Acel, temp, and uinterf, which have N rows and K columns
    prompt = ' Would you like to refine this further? (yes == 1) \n [Remember to rename SSNEW] ';
    redo = input(prompt);
end

global N K Gamma L P0 P0t    

% These are column vectors
% We import the steady state, with constant P0 
% Options for filenames: Regular parameter values with
    
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
data  = csvread('TextFiles/SSDataPeBi1Gamma1K300Qexp.csv');
A0    = data(2*K+1:3*K)'; 
u0bar = data(8*K+2:10*K+2)';
th0   =  data(4*K+1:6*K)';    
lam0  = data(10*K+3); 

% We define the sinusoidal Pin, with respect to the constant P0 used
% for the steady state
%% FOR SINUSOIDAL P 
tau      = linspace(0,2*pi,N)'; 
Pin      = P0t(tau/omega);
Pinbar   = mean(P0t(tau/omega));
Pintilde = Pin- Pinbar;
Pintau   = Pintilde.*tau;
% I don't need the analytical part but sometimes we use it for comparison
% if things look too weird
PinAn      = @(tau) P0+DeltaP*sin(tau); 
PintauAn   = @(tau) DeltaP*sin(tau).*tau; 
PintildeAn = @(tau) DeltaP*sin(tau);
PinbarAn   = @(tau) P0; 
% top , from 0 to lam0

xvec  = linspace(0,1,K);
xtop  = lam0.*xvec; 
dx    = (xtop(end)-xtop(1))/(K-1);
%bottom, from lam0 to 1
xbot  = lam0 + (L-lam0).*xvec; % for the integrals, this is okay, as lam0 is a constant.
% Integral from lam0 to x, where x goes from lam0 to 1 .The negative of
% this is the integral we want, i.e. integral from x to lam0 where x goes
% from lam0 to 1

%fac2 = derivative(A0.*Cint,dx); 
%tau = omegaP.*t;

phitau2  = zeros(N,1);
phitauAn = phitau2;

%tau = linspace(0,2*pi,N)';
dt = (1)/(N-1);
% analytic expression for Pintau
term1An = (1/(2*pi))*integral(PintauAn,0,2*pi);

%vector expression for Pintau
term1  = (1/(2*pi))*trapz(Pintau).*((2*pi-0)*dt); 

for i = 1:N
    phitauAn(i) = term1An + integral(PintildeAn,0,tau(i)); 
    phitau2(i)  = term1 + trapz(Pintilde(1:i))*(tau(i)-0)*(1/(i)); 
end
phitau= term1 + cumtrapz(Pintilde)*(2*pi-0).*dt; 
check = cumtrapz(Pintilde)*(2*pi-0).*dt; 

% We calculate all the ingredients for A1tilde
fth0    = flip(th0(1:K));
fA0     = flip(A0(1:K)); 
Cint    = flip(-cumtrapz(1./(exp(-Gamma*fth0).*fA0)).*(-1)*dx);  
fac1    = smooth(derivative(A0(1:K).*Cint,dx));
fac2    = smooth(derivative(A0(1:K),dx).*Cint + A0(1:K).*derivative(Cint,dx)); 
% Equivalent expressions for A1tilde (comment one of them)
%A1tilde = -(1/3)*phitau*fac2';
A1tilde = (1/3)*phitau*(exp(Gamma.*th0(1:K))-derivative(A0(1:K),dx).*Cint);
% We calculate t1tilde
dth0dx = smooth(derivative(th0(1:K),dx))';
t1tilde = -(1/3).*phitau*Cint.*(dth0dx);

% We use an anlytical expression for dAdx at lam0
A0primelam0 = (St.*trapz([A0(1:K), 1]).*(1).*dx + Pinbar)./(3*exp(-Gamma*th0(K))); 

% We calculate lam1tilde
lam1tilde = -phitau*(1./(3*A0primelam0*exp(-Gamma*th0(K))));

% We calculate u0tilde 
u0tilde = (Pintilde./3)*Cint(1:K);

% In order to plot the data, we have to interpolate. We define the matrices
% that allow  us to interpolate. 
tval = t(end)-2*pi/omega; 
% We shift so that everything is from 0 to 2 pi
indx = find(t>=tval,1);
indx = indx -1; 

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

xmatrix    = [lam*xcel',lam + (L-lam)*xcel'];
xmatrixint =  [lam*xint',lam + (L-lam)*xint(2:end)'];
xmat2      = lam*xcel';
xmat1int   = lam*xint';


Anum     = zeros(N,K);
thnum    = zeros(N,K);
u0barnum = zeros(N,K); 
for i = 1:N
    Anum(i,:)  = interp1(xmatrix(i,:),Acel(i,:),linspace(0,L,K),'pchip'); 
    thnum(i,:) = interp1(xmatrix(i,:),temp(i,:),linspace(0,L,K),'spline');
end

xforu = linspace(0,1,K+1)*lam0;
for i = 1:N
     u0barnum(i,:) = interp1(xmatrixint(i,:),uinterf(i,:),linspace(0,L,K),'spline'); 
end

% We compute the averages

omegat     = omega*t(indx:end);
AAvg       = mean(Anum(indx:end,:));
thAvg      = mean(thnum(indx:end,:));
u0Avg      = mean(u0barnum(indx:end,:));

% We subtract the averages

u0tildenum   = u0barnum - ones(N,1)*u0Avg; 

A1tildenum   = omega*(Anum-ones(N,1)*AAvg);

t1tildenum   = omega*(thnum-ones(N,1)*thAvg);  

lam1tildenum = omega*(lam-mean(lam(indx:end))); 


% We run the plotting code
run('PlottingFiles/ContoursLargeOmega1')


