% DATE:     2020 
% DESCR:    TimeDependentWithdPdtNoEvent
%           Main code for time dependent problem with CONTROL SYSTEM 1. Uses
%           ParametersDefinition, coupledPdeWithdPdt, and usolution to solve pdes.
%           Instead of using an event function, we just use a for loop to
%           separate the addition into days. This code is slower, but it
%           works and it is an alternative for anyone that doesn't want to
%           use inbuilt Events, although I highly recommend Events. 
%           Computation is stopped when P becomes smaller than 0.01. 
%           This is tol_P. Then new cylinders are added, according to some 
%           constant value or pattern, and the computation is restarted. 
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



%% COMPUTE STEADY STATE
% Define parameters
global N K T D uf P0 d
% We define all of the parameters in an external routine for clarity 
ParametersDefinition

%P0 = c1*D*St
% set to 1 if we want to compare with steady state
st = 1; 

if st == 1 
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

        A0      = (1- D).*linspace(0,1,K+1)'+D; 
        A0      = (A0(1:end-1)+A0(2:end))/2; 
        th0     = zeros(K,1); 
        phi0    = zeros(K,1);
        lam0    = 0.7;
        A0lamt  = lam0.*A0;
    case 'steady'
        A0     = data(2*K+1:3*K); 

        th0    =  data(4*K+1:5*K);
        phi0   = data(5*K+1:6*K); 
    
        lam0   = data(10*K+3); 
        A0lamt = A0.*lam0;
        P0     = data(10*K+4);
end


%Initialise the variables
Alamt   = []; 
th      = [];
phi     = [];
P       = [];
lam     = []; 

%% ODE integration 
te     = 0;
j      = 0;
tvec   = [];
cyladd = [];
first  = [];

prompt = ' Do you want to use a pattern addition ? (yes == 1)';
m      = input(prompt);


% if m == 1
%     fac = zeros(floor(T/d),1);
%     for i=1:floor(T/d)
%         val = (i) - 7*floor((i-1)/7);
%         if (val<=4)
%             fac(i) = 3;
%         elseif val== 5  
%             fac(i) = 16;
%         else 
%             fac(i) = 0; 
%         end
%     end
% else 
%     fac = (1/D).*ones(length(fac),1);
% end

tic
        y0(1:K) = A0lamt;
        y0(1+K:2*K) = A0.*th0;
        y0(2*K+1:3*K) = phi0;	
        y0(3*K+1) = P0; 
        y0(3*K+2) = lam0;  
      
tstart = 0;
        
for j=1:floor(T/d)
    disp(tstart)
    tend = j*d;
    tspan = [tstart, tend];
    first   = [first; tstart]; 
    %options = odeset('RelTol',1.0e-3,'AbsTol',1.0e-2);
    [t,y] = ode15s(@coupledPdeWithdPdt,tspan,y0);    
    Alamt   = [Alamt; y(2:end,1:K)]; % This is A from X=0 to X=1 (this is, 0<x<lambda)

    phi =  [phi; y(2:end,2*K+1:3*K)];

    P   = [ P; y(2:end,3*K+1)]; 

    lam = [lam; y(2:end,3*K+2)]; 

    Anow = y(2:end,1:K)./y(2:end,3*K+2);
    th  = [th; y(2:end,K+1:2*K)./Anow];

    tvec = [tvec; t(2:end)];

    y0(1:K) = y(end,1:K);
    y0(1+K:2*K) = y(end,K+1:2*K) ;
    y0(2*K+1:3*K) = y(end, 2*K+1:3*K);
    if m == 1
             val = j- 5*floor((j-1)/5);
             if val == 5  
                 fac = 3*c1;
             else 
                 fac = 1*c1;
             end
        else
        fac =  0.2800; 
    end
        

    y0(3*K+1) = y(end,3*K+1)+D*St*(fac); 
    cyladd        = [cyladd; fac];
    %y0(3*K+1) = ye(3*K+1)+D*St*(0.1286); 
    y0(3*K+2) = y(end,3*K+2);  

    tstart = tend; 

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
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i),1,P(i));   
end


% We add the Dirichlet boundary conditions 
thtop = [ zeros(N,1), ...
       th         ];

Acel = [ A, ones(N,K)];														% A (cell values)
Aint = ([D*ones(N,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)

uinterf  = [ u , ...
					uf.*ones(N,K+1) ];
temp = [th, phi];		% complete temperature profile (theta and phi)


xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
    
    
if st==0
    xvector1 = [xcel*lam(end);lam(end) + xcel*(L-lam(end))];
    xvector2 = [xint*lam(end);lam(end) + xint(2:end)*(L-lam(end))];
    SS = [xvector1; Acel(end,:)'; temp(end,:)'; xvector2; uinterf(end,:)'; lam(end); P(end)]; 
    csvwrite('SSData.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end

% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms
% set to 1 if we want to save data in csv file 
dat = 0; 
sav = 0;
P0tval = 0; 

%PlottingTimesteps
run('PlottingFiles/ContoursdPdt')
