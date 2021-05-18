% DATE:     2020 
% DESCR:    TimeDependentuftBimodalWithdPdt
%           Main code for time dependent problem with CONTROL SYSTEM 2 plus
%           a bimodal velocity function. This code takes a very long time
%           and is most likely not the right thing to do. It uses
%           ParametersDefinition, coupledPdeWithdPdt, and usolution to solve pdes.
%           New cylinders are added, according to some constant value or
%           pattern, and the computation is restarted. We fix the times at
%           which we add cylinders, but are chosen such that the
%           consumption rate is in balance with the addition rate. 
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
%           PlottingFiles/ContoursuftBimodal: Code where we plot and save the data for
%           the contours to illustrate the results computed here. 

% Clear previous files
close all 
clear all


%% COMPUTE STEADY STATE
% Define parameters
global N K T D P0t P0 uft

% We define all of the parameters in an external routine for clarity 
ParametersDefinition


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
    data = csvread('TextFiles/SSK300pattern2Qp7.csv');
end

P0  = 0.0412; % average P value for the Pattern case, CS 2a 
P0t = @(t) P0;

% Initial conditions

% incon can be steady to check return to steady, or simple, which is just
% linear A, th =0  everywhere 

incon = 'steady';

switch incon
    case 'simple'

        %A0 = (1- 1e-5 -D).*linspace(0,1,K)'+D; 
        
        A0   = (1- D).*linspace(0,1,K+1)'+D; 
        A0   = (A0(1:end-1)+A0(2:end))/2; 
        th0  = zeros(K,1); 
        phi0 = zeros(K,1);
        lam0 = 0.7;
        
    case 'steady'
        A0   = data(2*K+1:3*K); 

        th0  =  data(4*K+1:5*K);
        phi0 = data(5*K+1:6*K); 
    
        lam0 = data(10*K+3); 
end


%Initialise the variables
Alamt  = []; 
th     = [];
phi    = [];
P      = [];
lam    = []; 
u      = []; 
uftvec = [];
cyladd = [];

%% ODE integration variables 
te      = 0;
j       = 0;
tvec    = [];
first   = 0;
options = odeset('RelTol',1.0e-4,'AbsTol',1.0e-4);
onecyl  = D*St*c1;

y0(1:K)       = A0.*lam0;
y0(1+K:2*K)   = A0.*th0;
y0(2*K+1:3*K) = phi0;
y0(3*K+1)     = P0;
y0(3*K+2)     = lam0;  

% This code operates with a for loop which is in charge of
% the velocity changes. During each iteration of this group, we check the 
% level of P and add cylinders when necessary. The amount added corresponds
% to either a constant or a set pattern, which is already determined and 
% chosen such that it matches the consumption rate. 

%Independent variable for ODE integration 
tstart = 0;
%d = d/(32); 
dayfactor = 4;
d2 = d/dayfactor;
j2=0;

tic

T = floor(T/d).*d; % this enables us to finish the computation at a 
                   % specific day, rather than an awkward day fraction 
for j=1:floor(T/d2)
    disp(tstart)
    tend = j*d2;
%% ODE integration 
    % Changing within a day - unevenly distributed
    dayfrac = 3600*uc/(dayfactor*Ldim);
 
    uft = @(t) (24).*(t<=(dayfrac+tstart));

    tspan = [tstart tend] ; 
    
    [t,y] = ode15s(@coupledPdeWithdPdtuft,tspan,y0,options); 

    Alamt   = [Alamt; y(:,1:K)]; % This is A from X=0 to X=1 (this is, 0<x<lambda)

    phi =  [phi; y(:,2*K+1:3*K)];

    P = [P; y(:,3*K+1)]; 
    
    lam = [lam; y(:,3*K+2)]; 

    Anow = y(:,1:K)./y(:,3*K+2);
    th  = [th; y(:,K+1:2*K)./Anow];
    
    tvec = [tvec; t];
    
    y0(1:K) = y(end,1:K);
    y0(1+K:2*K) = y(end,K+1:2*K) ;
    y0(2*K+1:3*K) = y(end, 2*K+1:3*K);

    y0(3*K+2) = y(end,3*K+2);
    y0(3*K+1) = y(end,3*K+1);
    
    if mod(j,dayfactor)==0
        j2 = j2+1;
       Pend = y(end,3*K+1);
         
        val = j2- 7*floor((j2-1)/7);
            if val== 5   
              Pmax = uc*(24*60*60)*3/Ldim+ 0.01;
               %Pmax = 3*0.3456/Ldim + 0.01; 
               %Pmax = 3*0.6912/Ldim + 0.01;
              % Pmax = 3*0.5/Ldim + 0.01;
         
            elseif val==6 || val == 7
               Pmax = 0.01; 
            else 
               
               Pmax = uc*(24*60*60)/Ldim + 0.01 ;
               %Pmax = 0.3456/Ldim + 0.01;
              % Pmax = 0.6912/Ldim + 0.01;
               %Pmax = 0.5/Ldim +0.01; 
%             end
            end
             k=0;
       
       while (Pend<Pmax)
            % we can choose to add continous amounts of
            %paste, or in discrete cylinders
            
            %Pend = Pend+dPend;
            Pend = Pend+onecyl;
            k     = k+1;
        end
        %cyladd = [cyladd; k*dPend];
        cyladd    = [cyladd; k*onecyl]; 
        y0(3*K+1) = Pend;
    end

    tstart = tend; 

%% We calculate u with the solution for A, th and lam
    unew   = zeros(size(Anow));
    Ntemp  = length(t);
    for i=1:Ntemp
        unew(i,:) = usolutionuft(Anow(i,:)',(y(i,K+1:2*K)./Anow(i,:))',y(i,3*K+2),1,y(i,3*K+1),uft(t(i)));   
    end
    u      = [u; unew]; 
    %uftvec = [uftvec; uft(t)*ones(length(t),1)];
    uftvec = [uftvec; uft(t)];
end

toc 

N = length(tvec); 
t = tvec; 
% We add the Dirichlet boundary conditions 
thtop    = [ zeros(N,1), ...
           th         ];
A        = Alamt./lam; 
Acel     = [ A, ones(N,K)];														% A (cell values)
Aint     = ([2*D*ones(N,1) - A(:,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)
uinterf  = [ u , ...
					uftvec.*ones(N,K+1) ];
temp     = [th, phi];		% complete temperature profile (theta and phi)


xint     = linspace(0,1,K+1)';
xcel     = linspace(xint(2)/2,1-xint(2)/2,K)';
   

indx = length(t);
if st==0
    xvector1 = [xcel*lam(indx);lam(indx) + xcel*(L-lam(indx))];
    xvector2 = [xint*lam(indx);lam(indx) + xint(2:end)*(L-lam(indx))];
    SS = [xvector1; Acel(indx,:)'; temp(indx,:)'; xvector2; uinterf(indx,:)'; lam(indx); P0]; 
    csvwrite('SSData.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end

% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms
% set to 1 if we want to save data in csv file 
dat  = 0; 
sav  = 0;

uftval = 1 ; 
P0tval = 1; 

run('PlottingFiles/ContoursuftBimodal')
