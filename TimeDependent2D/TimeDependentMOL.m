% Clear previous files
close all 
clear all
clc


%% COMPUTE STEADY STATE
% Define parameters

% We define all of the parameters in an external routine for clarity 
ParametersDefinition

global N K T D uf P0t P0

options = odeset('RelTol',1.0e-03,'AbsTol',1.0e-06);

% set to 1 if we want to compare with steady state
st = 0; 

if st == 1 
    % Change filename to match what we want to import 
    data = csvread('SSG96K300.csv');
end

P0t = @(t)P0; 
%P0t = @(t) P0 + P0*sin(2*pi*t); % base case 
% Initial conditions

% incon can be steady to check return to steady, or simple, which is just
% linear A, th =0  everywhere 

incon = 'simple'; 

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

y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = phi0;	%must flip since th0bot is stored in reverse order
y0(3*K+1) = lam0;  

% Independent variable for ODE integration 
tout = linspace(0,T,N);

%% ODE integration 

tic
[t,y] = ode15s(@coupledPde,tout,y0); 
toc

A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

phi   = y(:,2*K+1:3*K);

lam = y(:,3*K+1); 

% Save the solution from here and then import into steady state to see if
% it actually converges to a steady state

%% We calculate u with the solution for A, th and lam
u = zeros(size(A));
for i=1:N
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i),1,P0t(tout(i)));   
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
    
    
if st==0
    xvector1 = [xcel*lam(end);lam(end) + xcel*(L-lam(end))];
    xvector2 = [xint*lam(end);lam(end) + xint(2:end)*(L-lam(end))];
    SS = [xvector1; Acel(end,:)'; temp(end,:)'; xvector2; uint(end,:)'; lam(end)]; 
    csvwrite('SSData.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end

% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms
% set to 1 if we want to save data in csv file 
dat = 0; 

%PlottingTimesteps
%PlottingContours
