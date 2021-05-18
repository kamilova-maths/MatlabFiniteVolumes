close all 
clear all

global N K T D uf P0 P0t
% We define all of the parameters in an external routine for clarity 
ParametersDefinition

%P0 = c1*D*St
% set to 1 if we want to compare with steady state
st = 1; 

if st == 1 
    % Change filename to matcth what we want to import 
    data = csvread('All1Omega100n105.csv');
    %    data = csvread('All1Gamma10.csv');
end


incon = 'steady'; 

switch incon
    case 'simple'

        A0 = (1- D).*linspace(0,1,K+1)'+D; 
        A0 = (A0(1:end-1)+A0(2:end))/2; 
        th0 = zeros(K,1); 
        phi0 = zeros(K,1);
        lam0 = 0.7;
        A0lamt = lam0.*A0;
    case 'steady'
        A0 = data(2*K+1:3*K); 

        th0 =  data(4*K+1:5*K);
        phi0 = data(5*K+1:6*K); 
    
        lam0 = data(10*K+3); 
        A0lamt = A0.*lam0;
        %P0   = data(10*K+4);
end

%% START FROM HERE WHEN YOU HAVE ALREADY CALCULATED STEADY STATE 

% Initial conditions 
% 
% A0 = A0steady; 
% th0 = th0steady; 
% th0bot = phi0steady; 

%Q = Bi; % can my code cope with that? Even if the steady state code can't. 
%[Probably not but sort of worth a chance]
y0(1:K) = A0lamt;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = phi0;		
y0(3*K+1) = lam0; 

%P0t = @(t)P0*cos(2*pi*t/0.1234);
%P0 = 1; 
omega = 100;
DeltaP = 0.8/P0; 
n = 105;
T= 2*pi*n/omega;
%T = 10; 

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

Acel = [ A, ones(N,K)];														% A (cell values)
Aint = ([D*ones(N,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)

uint  = [ u , ...
					uf.*ones(N,K+1) ];
temp = [th, phi];		% complete temperature profile (theta and phi)

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
    
st=0;    
if st==0
    xvector1 = [xcel*lam(end);lam(end) + xcel*(L-lam(end))];
    xvector2 = [xint*lam(end);lam(end) + xint(2:end)*(L-lam(end))];
    SS = [xvector1; Acel(end,:)'; temp(end,:)'; xvector2; uint(end,:)'; lam(end-1); P0t(t(end))]; 
    csvwrite('NewFile.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end

% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms
% set to 1 if we want to save data in csv file 
dat = 0; 
sav = 0; % indicator for saving data
P0tval = 1;
uftval = 0;
%PlottingTimesteps
%PlottingContoursOmegat