% Clear previous files
close all 
clear all
clc


%% COMPUTE STEADY STATE
% Define parameters

% We define all of the parameters in an external routine for clarity 
ParametersDefinition

global N K T D uf dCdt P0 d counter


% set to 1 if we want to compare with steady state
st = 1; 

if st == 1 
    % Change filename to match what we want to import 
    data = csvread('SSG125K300Qorg.csv');
end

%P0t = @(t)P0; 
%P0t = @(t) P0 + P0*sin(2*pi*t); % base case 
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

y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = phi0;	%must flip since th0bot is stored in reverse order
y0(3*K+1) = P0;
y0(3*K+2) = lam0;  

% Independent variable for ODE integration 
tout = linspace(0,T,N);

% We define non-dimensional day d
d = 86400*uc/Ldim; 
counter = 0; 
% We define the addition of cylinders

%dCdt = @(t) 0; 
% term1 =@(t)0;
% term2 =@(t)0;
% term3= @(t)0;
% for j =1:floor(T/d)
%    term1 = @(t)term1 + dirac(t-j*d); 
% end
% 
% for j =0:floor(T/(3*d))
%    term2 = @(t)term2 + dirac(t-(2+3*j)*d); 
% end
% 
% for j = 1:floor(T/(3*d))
%     term3= @(t)term3 + dirac(t-3*j*d); 
% end


freq = 1/d;
%dCdt = @(t)term1(t) + term2(t) + term3(t); 

%dCdt =@(t)dirac(t-d)+dirac(t-2*d)+dirac(t-2*d); 
% vec=ones(size(tout));
% fac = zeros(floor(T/d),1); 
% for j=1:floor(T/d)
%     fac(j) = j - 3*floor((j-1)/3);
%     vec(j*N/floor(T/d)) = 0 ;
% end
% dCdtv = dirac(vec); 
% % Replace peaks by factors
% idx = C == Inf; % find Inf
% dCdtv(idx) = fac;   % You have to be a function! [>.<]
% WEIRD ATTEMPT - PLEASE CHANGE AFTERWARDS (IF CODE RUNS)
C1=@(t) 1.5*sawtooth(2*pi*t*8)+1.5; 
C2=@(t) sawtooth(2*pi*t*4)+1; 
C1=@(t) sawtooth(2*pi*(t-1/8)*8) + 1; 
%dCdt = @(t)C1(t) -C2(t) +1; 
%dCdt = @(t)sawtooth(2*pi*t*8) + 2;
%dCdt = 0; 
dCdt = @(t)0.5*sawtooth(2*pi*t*freq) + 0.5; 
%u0t = 4; 
%Cdt = derivative(C,T/N);
%% ODE integration 

%for i=1:floor(T/d)
%options = odeset('RelTol',1.0e-03,'AbsTol',1.0e-06,'Events',@EventFunction);
options = odeset('RelTol',1.0e-03,'AbsTol',1.0e-06);

tic
%[t,y,te,ye,ie] = ode15s(@coupledPdeWithdPdt,tout,y0,options); 
[t,y] = ode15s(@coupledPdeWithdPdt,tout,y0,options); 
toc

A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

phi   = y(:,2*K+1:3*K);

P = y(:,3*K+1); 

lam = y(:,3*K+2); 

%end

% Save the solution from here and then import into steady state to see if
% it actually converges to a steady state

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

PlottingTimesteps
%PlottingContours
