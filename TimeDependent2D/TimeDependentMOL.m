% DATE:     2020 
% DESCR:    TimeDependentMOL
%           Main code for time dependent problem. Uses
%           ParametersDefinition, coupledPde, and usolution to solve pdes. 
%
% INPUT: 
%           No input variables
%          
% OUTPUT:   Main outcomes: 
%           Acel: The N x 2*K matrix storing the cell values for A 
%           uinterf: The N x (2*K +1) matrix storing all interface values for
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
%           PlottingFiles/Timesteps: Plots solutions at particular timesteps.
%           Also has the option to compare with steady state and save the
%           .csv files.


% Clear previous files
close all 
clear all


%% COMPUTE STEADY STATE
% Define parameters
global P0t
ParametersDefinition
% We define all of the parameters in an external routine for clarity 

options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-8);

% set to 1 if we want to compare with steady state
st = 0; 

if st == 1 
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
    % PeBiGamma1Omega100Qexp
    % PeBi27Omega100Qexp
    % PeBi27Omega300Qexp
    data = csvread('TextFiles/SSDataP01.csv');
end

% We set the value for P0t, a time dependent version of imposed P0. If you
% would like to compute the constant P0 value, set the function to be equal
% to P0 for all t. 
prompt = ' Would you like constant P0 (input == 0) or sinusoidal P0 (input == 1) ? ';
p = input(prompt);

% We can either use simple conditions (which satisfy the boundary
% conditions) to evolve our solution to a steady state, or from a nearby
% steady state, to refine a current understanding of steady state. 

if p == 0
    prompt = ' What value of P0 would you like to set ? ';
    P0  = input(prompt);
    P0t =@(t) P0;
    incon  = 'simple';
elseif p == 1
    prompt = ' What value of P0 would you like to set ? ';
    P0     = input(prompt);
    prompt = ' Amplitude for P0 (as fraction of P0 )';
    DeltaP = input(prompt)/P0;
    prompt = ' Frequency for P0 ' ;
    omega  = input(prompt);
    prompt = ' How many periods shall we run ? ';
    n      = input(prompt);
    T      = 2*pi*n/omega;
    P0t    = @(t) P0 + DeltaP.*sin(omega*t); % base case 
    incon  = 'steady';
end

%incon = 'simple';

switch incon
    case 'simple'
        th0  = zeros(K,1); 
        phi0 = zeros(K,1);
        lam0 = 0.7;
        A0   = (1- D).*linspace(0,lam0,K)'./lam0+D; 
    case 'steady'
        A0   = data(2*K+1:3*K); 

        th0  =  data(4*K+1:5*K);
        phi0 = data(5*K+1:6*K); 
    
        lam0 = data(10*K+3); 
end

% We store the initial values. Note that we solve for A times lambda and 
% A times theta, splitting up theta between top and bottom. 

y0(1:K)       = A0.*lam0;
y0(1+K:2*K)   = A0.*th0;
y0(2*K+1:3*K) = phi0;	
y0(3*K+1)     = lam0;  

% Independent variable for ODE integration 
prompt = 'tout (l=1) or tspan (otherwise)? ';
l = input(prompt);

tic
if l == 1
    tout  = linspace(0,T,N);
    [t,y] = ode15s(@coupledPde,tout,y0); 
else
    tspan = [0 T] ; 
    [t,y] = ode15s(@coupledPde,tspan,y0); 
end

toc

%% ODE integration 

N     = length(t); 

Alam  = y(:,1:K); % This is A times lambda from X=0 to X=1 (this is, 0<x<lambda)

phi   = y(:,2*K+1:3*K);

lam   = y(:,3*K+1); 

A     = Alam./lam;

th    = y(:,K+1:2*K)./A;

%% We calculate u with the solution for A, th and lam

u = zeros(size(A));
for i=1:N
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i),1,P0t(t(i)));   
end

% We add the Dirichlet boundary conditions 
thtop = [ zeros(N,1), ...
       th         ];

Acel = [ A, ones(N,K)];									% A (cell values)
Aint = ([2*D*ones(N,1) - A(:,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)
uinterf = [ u , ...
					uf.*ones(N,K+1) ];
temp = [th, phi];		% complete temperature profile (theta and phi)

% Interface values of x, for 0<x<lambda
xint = linspace(0,1,K+1)';

% cell vaues for x, for 0<x<lambda
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';


% set to 1 if we want to save data in csv file 
dat = 1; 

% set to 1 if we want to save the png figures
sav= 0;

% set P0tval == 0 if P is a solution to the pde instead of imposed. Then the
% result for P will be plotted with respect to t, together with lambda. 
% set P0tval == 1 if P0t was imposed and you would like to plot it
% set P0tval = other value if you do not want to print P0t at all. 

P0tval = 2; 

% set uftval == 1 if we also solve for uf(t) and would like to plot it. Set
% to any other value otherwise

uftval =0;
prompt = 'Do you want to plot stuff ? (yes == 1) \n ';
j = input(prompt);
if j ==1 
    %run('PlottingFiles/ContoursOrg')
    run('PlottingFiles/Timesteps')
end

% Storing a new steady state text file, or update an existing one
prompt = ' Would you like to store this solution ? (yes == 1 ) [Stored as SSNEW] \n ';

st = input(prompt);   
if st==1
    xvector1 = [xcel*lam(end);lam(end) + xcel*(L-lam(end))];
    xvector2 = [xint*lam(end);lam(end) + xint(2:end)*(L-lam(end))];
    SS       = [xvector1; Acel(end,:)'; temp(end,:)'; xvector2; uinterf(end,:)'; lam(end); P0]; 
    csvwrite('TextFiles/SSNEW.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end

