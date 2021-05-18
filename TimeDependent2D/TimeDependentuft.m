% DATE:     2020 
% DESCR:    TimeDependentuft
%           Main code for time dependent problem imposing a time dependent
%           uft(t) . Uses ParametersDefinition, coupledPdeuft, and 
%           usolution to solve pdes. 
%
% INPUT: 
%           No input variables
%          
% OUTPUT:   Main outcomes: 
%           Acel: The N x 2*K matrix storing the cell values for A 
%           uint: The N x (2*K +1) matrix storing all interface values for
%           u
%           temp: The N x (2*K) matrix that stores cell values for
%           temperature
% ADDITIONAL COMMENTS: 
%
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           coupledPde: This is where the pdes are, that will be then
%           solved with method of lines with ode15s
%           PlottingFiles/ContoursOrg: Code where we plot and save the data for
%           the contours to illustrate the results computed here. 
%          

% Clear previous files
close all 
clear all


%% COMPUTE STEADY STATE
% Define parameters

global P0t uft
% We define all of the parameters in an external routine for clarity 
ParametersDefinition


options = odeset('RelTol',1.0e-4,'AbsTol',1.0e-4);

% set to 1 if we want to compare with steady state
st = 1; 

if st == 1 
    % Change filename to match what we want to import 
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
    % SSDatauftomega2Deltaup5.csv: P0 = 1, K = 300, uft with omega 2 and
    % Deltau 0.5
    data = csvread('TextFiles/SSDatauftomega2Deltaup5.csv');
end

P0t =@(t) P0;

% We define uft as a sinusoidal function 
omega = 2; 
Deltau = 0.5;
% to make the extraction of data easier at the end, T will be in terms of
% this omega
n = 5; % how many periods do we want to run this for
T= 2*pi*n/omega; 
uft = @(t) uf + Deltau*sin(omega*t);
% Initial conditions

% incon can be steady to check return to steady, or simple, which is just
% linear A, th =0  everywhere 

incon = 'steady';

switch incon
    case 'simple'
        A0 = (1- D).*linspace(0,1,K+1)'+D; 
        A0 = (A0(1:end-1)+A0(2:end))/2; 
        th0 = zeros(K,1); 
        phi0 = zeros(K,1);
        lam0 = 0.7;
        
    case 'steady'
        A0 = data(2*K+1:3*K); 

        th0 =  data(4*K+1:5*K);
        phi0 = data(5*K+1:6*K); 
    
        lam0 = data(10*K+3); 
end

y0(1:K) = A0.*lam0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = phi0;
y0(3*K+1) = lam0;  

% Independent variable for ODE integration 

%% ODE integration 
l = 0;
tic
if l == 1
    tout  = linspace(0,T,N);
    [t,y] = ode15s(@coupledPdeuft,tout,y0); 
else
    tspan = [0 T] ; 
    [t,y] = ode15s(@coupledPdeuft,tspan,y0); 
    N = length(t); 
end

toc

Alam  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

phi   = y(:,2*K+1:3*K);

lam = y(:,3*K+1); 

A = Alam./lam;

th = y(:,K+1:2*K)./A;


%% We calculate u with the solution for A, th and lam
u = zeros(size(A));
for i=1:N
    u(i,:) = usolutionuft(A(i,:)',th(i,:)',lam(i),1,P0t(t(i)),uft(t(i)));   
end


% We add the Dirichlet boundary conditions 
thtop = [ zeros(N,1), ...
       th         ];

Acel = [ A, ones(N,K)];		% A (cell values)
Aint = ([2*D*ones(N,1) - A(:,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)
uint  = [ u , ...
					uft(t).*ones(N,K+1) ];
temp = [th, phi];		% complete temperature profile (theta and phi)


xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
    
indx = N;
if st==0
    xvector1 = [xcel*lam(indx);lam(indx) + xcel*(L-lam(indx))];
    xvector2 = [xint*lam(indx);lam(indx) + xint(2:end)*(L-lam(indx))];
    SS = [xvector1; Acel(indx,:)'; temp(indx,:)'; xvector2; uint(indx,:)'; lam(indx); P0]; 
    csvwrite('TextFiles/SSDatauftomega2Deltau5.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end

% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms
% set to 1 if we want to save data in csv file 
dat = 0; 

P0tval = 2; 
uftval = 1; 

prompt = 'Do you want to plot stuff ? (yes == 1) \n ';
j = input(prompt);
if j ==1 
    run('PlottingFiles/ContoursOrg')
    %run('PlottingFiles/Timesteps')
end
