% Clear previous files
close all 
clear all
clc


%% COMPUTE STEADY STATE
% Define parameters
global N K T D uf P0 c1 d
% We define all of the parameters in an external routine for clarity 
ParametersDefinition

%P0 = c1*D*St
% set to 1 if we want to compare with steady state
st = 1; 

if st == 1 
    % Change filename to match what we want to import 
    %data = csvread('SSK300avg2.csv');
    data = csvread('SSK300Isothermal2.csv');
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
        P0   = data(10*K+4);
end


% Independent variable for ODE integration 
tout = linspace(0,T,N);




% We define the addition of cylinders

%Initialise the variables
Alamt   = []; 
th  = [];
phi = [];
P   = [];
lam = []; 
first = [];

%% ODE integration 
te=0;
j=0;
tvec=[];

cyladd = []; 

tic

while (isempty(te)==0)
    
    if j==0
        y0(1:K) = A0lamt;
        y0(K+1) = P0; 
        y0(K+2) = lam0;  
        tspan = [0 T];
    else
             val = j- 7*floor((j-1)/7);
             if val== 5  
                 fac = 20*c1;
             else 
                 fac = 2*c1;
             end
        %fac =  0.4700; 
        cyladd = [cyladd; fac];
        y0(1:K) = ye(1:K);
        y0(K+1) = ye(K+1)+D*St*(fac); 
        y0(K+2) = ye(K+2);  

       tspan = [te T];
    end
j=j+1;

options = odeset('RelTol',1.0e-4,'AbsTol',1.0e-4,'Events',@EventFunctionIsothermal);

% te - column vector of the times at which events occurred
% ye - contains the solution value at each of the event times in t.e.
% ie - contains indices into the vector returned by the event function. The
% values indicate which event the solver detected
[t,y,te,ye,ie] = ode15s(@coupledPdeWithdPdtIsothermal,tspan,y0,options); 
%[t,y] = ode15s(@coupledPdeWithdPdt,tout,y0,options); 
te
first = [first; te]; 
%The 't' values are given at the rows 

Alamt   = [Alamt; y(:,1:K)]; % This is A from X=0 to X=1 (this is, 0<x<lambda)

P   = [ P; y(:,K+1)]; 
 
lam = [lam; y(:,K+2)]; 

%A = y(:,1:K)./y(:,K+2);

tvec = [tvec; t];

end
cylvec = [first, cyladd] ; 
csvwrite('DStCyl.csv', cylvec); 
A = Alamt./lam; 
avg = sum(cyladd)/length(first) 
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

Acel = [ A, ones(N,K)];														% A (cell values)
Aint = ([D*ones(N,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)

uint  = [ u , ...
					uf.*ones(N,K+1) ];
temp = zeros(N,2*K);

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

%indx = first(end);
indx = find(t==first(end-1),1);
indx=indx+1;
if st==0
    xvector1 = [xcel*lam(indx);lam(indx) + xcel*(L-lam(indx))];
    xvector2 = [xint*lam(indx);lam(indx) + xint(2:end)*(L-lam(indx))];
    SS = [xvector1; Acel(indx,:)'; temp(indx,:)'; xvector2; uint(indx,:)'; lam(indx); P(indx)]; 
    csvwrite('SSData.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end

% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms
% set to 1 if we want to save data in csv file 
dat = 0; 
sav = 0; % indicator for saving data
P0tval =0;
%PlottingTimesteps
PlottingContours