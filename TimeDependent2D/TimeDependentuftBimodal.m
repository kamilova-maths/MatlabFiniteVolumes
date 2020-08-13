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
    % Change filename to match what we want to import 
    data = csvread('SSK300uft2.csv');
end

P0 = 1; 
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
tstart = 0;
d = d/2;
for j=1:floor(T/d)
    tstart
    tend = j*d;
%% ODE integration 
    if mod(j,2)==0
        uft = @(t) 1.5;
    else
        uft = @(t) 0.5;
    end
    tic
    tspan = [tstart tend] ; 
    [t,y] = ode15s(@coupledPdeuft,tspan,y0); 
    toc

    Alamt   = [Alamt; y(:,1:K)]; % This is A from X=0 to X=1 (this is, 0<x<lambda)

    phi =  [phi; y(:,2*K+1:3*K)];

    lam = [lam; y(:,3*K+1)]; 

    Anow = y(:,1:K)./y(:,3*K+1);
    th  = [th; y(:,K+1:2*K)./Anow];
    
    tvec = [tvec; t];

    y0(1:K) = y(end,1:K);
    y0(1+K:2*K) = y(end,K+1:2*K) ;
    y0(2*K+1:3*K) = y(end, 2*K+1:3*K);

    y0(3*K+1) = y(end,3*K+1); 
  

    tstart = tend; 

%% We calculate u with the solution for A, th and lam
    unew = zeros(size(Anow));
    Ntemp = length(t);
    for i=1:Ntemp
        unew(i,:) = usolutionuft(Anow(i,:)',(y(i,K+1:2*K)./Anow(i,:))',y(i,3*K+1),1,P0t(t(i)),uft(t(i)));   
    end
    u = [u; unew]; 
    uftvec = [uftvec; uft(t)*ones(length(t),1)];
end

N = length(tvec); 
t = tvec; 
% We add the Dirichlet boundary conditions 
thtop = [ zeros(N,1), ...
       th         ];
A = Alamt./lam; 
Acel = [ A, ones(N,K)];														% A (cell values)
Aint = ([2*D*ones(N,1) - A(:,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)
uint  = [ u , ...
					uftvec.*ones(N,K+1) ];
temp = [th, phi];		% complete temperature profile (theta and phi)


xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
    
indx = find(t>pi,1); 
indx = indx +1; 
if st==0
    xvector1 = [xcel*lam(indx);lam(indx) + xcel*(L-lam(indx))];
    xvector2 = [xint*lam(indx);lam(indx) + xint(2:end)*(L-lam(indx))];
    SS = [xvector1; Acel(indx,:)'; temp(indx,:)'; xvector2; uint(indx,:)'; lam(indx); P0]; 
    csvwrite('SSData.csv', SS); 
    disp('Remember to change the name of the file at the end. Include Gamma and K')
end

% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms
% set to 1 if we want to save data in csv file 
dat = 0; 
sav= 0;

uftval = 1 ; 
% plot(t,u(:,1),t,1./Aint(:,1))
% hold on 
% plot(t,lam)
% 
% u(end,1)
%PlottingTimesteps
PlottingContoursuft
