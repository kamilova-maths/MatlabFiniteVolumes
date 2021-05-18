% DATE:     2020 
% DESCR:    TimeDependentuftBimodal
%           Main code for time dependent problem with a bimodal expression 
%           for uft. Here, uft is either 'on' or 'off', in such a way that
%           the daily average is 1. The velocity is turned 'on' for a
%           certain amount of time, a dayfrac, which is chosen such that
%           the average is still one when imposing the 'on' value as 24.
%           This was the simplest way for me to control the bimodal
%           variation, so, fixed magnitude of uft max, fixed average uft at
%           1, and we vary the dayfraction and the amount of 'drops'. 
%           Uses ParametersDefinition, coupledPde, and usolution to solve pdes. 
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
%           PlottingFiles/ContoursuftBimodal: Code where we plot and save the data for
%           the contours to illustrate the results computed here, bimodal uft. 
%           PlottingFiles/Timesteps: Plots solutions at particular timesteps.
%           Also has the option to compare with steady state and save the
%           .csv files.


% Clear previous files
close all 
clear all


%% COMPUTE STEADY STATE
% Define parameters
global N K T D P0t P0 uft

% We define all of the parameters in an external routine for clarity 
ParametersDefinition

options = odeset('RelTol',1.0e-4,'AbsTol',1.0e-4);

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
    data = csvread('TextFiles/SSK300P0p5.csv');
end

% This had been used for comparisons with the control systems, so one
% example is to use the average P value for the Pattern case, CS 2a. 
P0 = 0.0412; % average P value for the Pattern case, CS 2a 

P0t =@(t) P0;
%P0t = @(t) P0 + 0.7*P0*sin(2*pi*t); % base case 
% Initial conditions

% incon can be steady to check return to steady, or simple, which is just
% linear A, th =0  everywhere 

incon = 'steady';

switch incon
    case 'simple'

        %A0 = (1- 1e-5 -D).*linspace(0,1,K)'+D; 
        
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


%Initialise the variables
Alamt   = []; 
th  = [];
phi = [];
P   = [];
lam = []; 
u = []; 
uftvec = [];
%% ODE integration 
te=0;
j=0;
tvec=[];
first = 0;

y0(1:K) = A0.*lam0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = phi0;
y0(3*K+1) = lam0;  

% Independent variable for ODE integration 
% We modify the provided T so that it ends exacly at a integer number of
% day. This is in order to facilitate the plotting in ContoursuftBimodal 
T = floor(T/d).*d; 
tstart = 0;

% We have decided to split up the day in dayfactor amounts, and each amount
% will then have a section where uft is 'on' and a larger section where uft
% is 'off'. 

dayfactor = 4*3;
% We split up the days into 'dayfactor; chunks

d2 = d/dayfactor;
tic
for j=1:floor(T/d2)
    
    % We display tstart at each iteration as a way to track progress when
    % running the code 
    
    disp(tstart)
    tend = j*d2;
    % Changing within a day - unevenly distributed
    % dayfrac is chosen in such a way that for whatever dayfactor we
    % choose, we will allow the velocity to increase for a certain amount
    % of time ,corresponding to dayfrac, such that the average within a day
    % is still 1. We choose uft to be either 24 or zero, so the dayfrac is
    % such that when you average that out, it will be 1. 
    
    dayfrac = 3600*uc/(dayfactor*Ldim);
    %uft = @(t) 2.*(t>(tend-dayfrac));
    %uft = @(t) (540).*(t<=(dayfrac+tstart));
    % When I change uc instead of uf
    uft = @(t) (24).*(t<=(dayfrac+tstart));
    %uft = @(t) 2.*(t>=(tend-dayfrac)) + 0.9565.*(t<(tend-dayfrac));
    % A check
    %uft = @(t) 1;

    tspan = [tstart tend] ; 
    [t,y] = ode15s(@coupledPdeuft,tspan,y0);  

    Alamt   = [Alamt; y(2:end,1:K)]; % This is A from X=0 to X=1 (this is, 0<x<lambda)

    phi =  [phi; y(2:end,2*K+1:3*K)];

    lam = [lam; y(2:end,3*K+1)]; 

    Anow = y(2:end,1:K)./y(2:end,3*K+1);
    th  = [th; y(2:end,K+1:2*K)./Anow];
    
    tvec = [tvec; t(2:end)];
    % We use the last value as the initial value for the re-start
    y0(1:K) = y(end,1:K);
    y0(1+K:2*K) = y(end,K+1:2*K) ;
    y0(2*K+1:3*K) = y(end,2*K+1:3*K);

    y0(3*K+1) = y(end,3*K+1); 
  

    tstart = tend; 

%% We calculate u with the solution for A, th and lam
    unew = zeros(size(Anow));
    Ntemp = length(t(2:end));
    for i=1:Ntemp
        unew(i,:) = usolutionuft(Anow(i,:)',(y(i,K+1:2*K)./Anow(i,:))',y(i,3*K+1),1,P0t(t(i)),uft(t(i)));   
    end
    u = [u; unew]; 
    %uftvec = [uftvec; uft(t)*ones(length(t),1)];
    uftvec = [uftvec; uft(t(2:end))];
end
toc

N = length(tvec); 
t = tvec; 
% We add the Dirichlet boundary conditions 
thtop = [ zeros(N,1), ...
       th         ];
A = Alamt./lam; 
Acel = [ A, ones(N,K)];								% A (cell values)
Aint = ([2*D*ones(N,1) - A(:,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)
uinterf  = [ u , ...
					uftvec.*ones(N,K+1) ];
temp = [th, phi];		% complete temperature profile (theta and phi)


xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
    
% Save the data, if you need to use the.csv file of the last timestep
indx = length(t);
if st==0
    xvector1 = [xcel*lam(indx);lam(indx) + xcel*(L-lam(indx))];
    xvector2 = [xint*lam(indx);lam(indx) + xint(2:end)*(L-lam(indx))];
    SS = [xvector1; Acel(indx,:)'; temp(indx,:)'; xvector2; uinterf(indx,:)'; lam(indx); P0]; 
    csvwrite('SSData.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end


uftval = 1 ; 
P0tval = 2;

run('PlottingFiles/ContoursuftBimodal')
