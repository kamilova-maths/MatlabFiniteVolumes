% DATE:     2020 
% DESCR:    TimeDependentWithdPdtIsothermal
%           Main code for time dependent problem with CONTROL SYSTEM 1, for 
%           the isothermal case. Uses
%           ParametersDefinition, coupledPdeWithdPdtIsothermal, and usolution to solve pdes.
%           An Event function is used to stop computation, which checks for
%           when P becomes smaller than 0.01. This is tol_P. Then new
%           cylinders are added, according to some constant value or
%           pattern, and the computation is restarted. 
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
%           coupledPdeWithdPdt: This is where the pdes are, that will be then
%           solved with method of lines with ode15s. This includes the pde
%           for  P(t). 
%           PlottingFiles/Contours: Code where we plot and save the data for
%           the contours to illustrate the results computed here. 

% Clear previous files
close all 
clear all
clc


%% COMPUTE STEADY STATE
% Define parameters
global N K T D uf P0 c1 d St
% We define all of the parameters in an external routine for clarity 
ParametersDefinition
%St = 27;
%P0 = c1*D*St
% set to 1 if we want to compare with steady state
st = 1; 

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
    data = csvread('TextFiles/SSDataP01.csv');  
end


% We can either use simple conditions (which satisfy the boundary
% conditions) to evolve our solution to a steady state, or from a nearby
% steady state, to refine a current understanding of steady state. 


incon = 'steady'; 

switch incon
    case 'simple'

        A0     = (1- D).*linspace(0,1,K+1)'+D; 
        A0     = (A0(1:end-1)+A0(2:end))/2; 
        th0    = zeros(K,1); 
        phi0   = zeros(K,1);
        lam0   = 0.7;
        A0lamt = lam0.*A0;
    case 'steady'
        A0     = data(2*K+1:3*K); 

        th0    =  data(4*K+1:5*K);
        phi0   = data(5*K+1:6*K); 
    
        lam0   = data(10*K+3); 
        A0lamt = A0.*lam0;
        P0     = data(10*K+4);
end

% We define the addition of cylinders

%Initialise the variables
Alamt   = []; 
th      = [];
phi     = [];
P       = [];
lam     = []; 
first   = [];

%% ODE integration 

te     = 0;
j      = 0;
tvec   = [];

cyladd = []; 

options = odeset('RelTol',1.0e-4,'AbsTol',1.0e-4,'Events',@EventFunctionIsothermal);

% Select either pattern addition or constant addition 

prompt = 'Input 1 for pattern addition, 0 otherwise \n';
m      = input(prompt);


tic

while (isempty(te)==0)
    
    if j==0
        y0(1:K) = A0lamt;
        y0(K+1) = P0; 
        y0(K+2) = lam0;  
        tspan = [0 T];
    else
        if m ==1 
             val = j- 5*floor((j-1)/5);
             if val== 5  
                 fac = 3*c1;
             else 
                 fac = 1*c1;
             end
        else
             fac =  0.27826; % this is the average value (I think) 
        end
        cyladd   = [cyladd; fac];
        y0(1:K)  = ye(1:K);
        y0(K+1)  = ye(K+1)+D*St*(fac); 
        y0(K+2)  = ye(K+2);  

       tspan     = [te T];
    end
    
    j=j+1;


    % te - column vector of the times at which events occurred
    % ye - contains the solution value at each of the event times in t.e.
    % ie - contains indices into the vector returned by the event function. The
    % values indicate which event the solver detected
    [t,y,te,ye,ie] = ode15s(@coupledPdeWithdPdtIsothermal,tspan,y0,options); 
    %[t,y] = ode15s(@coupledPdeWithdPdt,tout,y0,options); 
    disp(te)
    first = [first; te]; 
    %The 't' values are given at the rows 

    Alamt   = [Alamt; y(2:end,1:K)]; % This is A from X=0 to X=1 (this is, 0<x<lambda)

    P   = [ P; y(2:end,K+1)]; 

    lam = [lam; y(2:end,K+2)]; 

    tvec = [tvec; t(2:end)];

end

A = Alamt./lam;
toc

% Save the solution from here and then import into steady state to see if
% it actually converges to a steady state
t = tvec; 
N = length(tvec); 
%% We calculate u with the solution for A, th and lam
u = zeros(size(A));
for i=1:N
    u(i,:) = usolution(A(i,:)',zeros(K,1),lam(i),1,P(i));   
end

Acel  = [ A, ones(N,K)];					% A (cell values)
Aint  = ([D*ones(N,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)

uinterf  = [ u , ...
		uf.*ones(N,K+1) ];
temp  = zeros(N,2*K);

xint  = linspace(0,1,K+1)';
xcel  = linspace(xint(2)/2,1-xint(2)/2,K)';

%indx = first(end);
indx = find(t==first(end-1),1);
indx=indx+1;
if st==0
    xvector1 = [xcel*lam(indx);lam(indx) + xcel*(L-lam(indx))];
    xvector2 = [xint*lam(indx);lam(indx) + xint(2:end)*(L-lam(indx))];
    SS = [xvector1; Acel(indx,:)'; temp(indx,:)'; xvector2; uinterf(indx,:)'; lam(indx); P(indx)]; 
    csvwrite('SSData.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end

P0tval = 0;
uftval = 0; 

run('PlottingFiles/ContoursdPdt')
%PlottingTimesteps
%PlottingContours
