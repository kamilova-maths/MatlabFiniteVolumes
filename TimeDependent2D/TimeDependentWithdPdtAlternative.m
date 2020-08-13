% Clear previous files
close all 
clear all
clc


%% COMPUTE STEADY STATE
% Define parameters
global N K T D uf P0 c1 d St Ldim
% We define all of the parameters in an external routine for clarity 
ParametersDefinition

%P0 = c1*D*St
% set to 1 if we want to compare with steady state
st = 1; 

if st == 1 
    % Change filename to match what we want to import 
    %data = csvread('SSApproach2PatternAfterPeak.csv');
    data = csvread('SSData.csv');
end

%P0 = D*St*0.1286;
%P0 = D*St*c1;


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
        A0lamt = lam0.*A0;
    case 'steady'
        A0 = data(2*K+1:3*K); 

        th0 =  data(4*K+1:5*K);
        phi0 = data(5*K+1:6*K); 
    
        lam0 = data(10*K+3); 
        A0lamt = A0.*lam0;
        P0   = data(10*K+4);
end

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
tstart = 0; 

% We define the addition of cylinders

cyladd = []; 
pattern = [1; 1; 1; 1; 3; 0; 0]; 
pattern = repmat(pattern,floor(floor((T+1)/d)/length(pattern)),1);
onecyl = D*St*c1; 
pattern = onecyl*pattern; 

%avgvalue = 0.1688;
%P0 = 4*onecyl;
N = 80;
tic 
dPend = 0.001; 

for j=1:floor(T/d)
    tend = j*d
    if j==1
        y0(1:K) = A0lamt;
        y0(1+K:2*K) = A0.*th0;
        y0(2*K+1:3*K) = phi0;	
        y0(3*K+1) = P0; 
        y0(3*K+2) = lam0;  
        %tspan = [0 tend];
        tout = linspace(tstart,tend,N);
    else
        Pend = y(end,3*K+1);
         
        val = j- 7*floor((j-1)/7);
            if val== 5   
              %Pmax = uc*(24*60*60)*3/Ldim+ 0.01;
               %Pmax = 3*0.3456/Ldim + 0.01; 
               %Pmax = 3*0.6912/Ldim + 0.01;
               Pmax = 3*0.5/Ldim + 0.01;
         
            elseif val==6 || val == 7
               Pmax = 0.01; 
            else 
               
               %Pmax = uc*(24*60*60)/Ldim + 0.01 ;
               %Pmax = 0.3456/Ldim + 0.01;
              % Pmax = 0.6912/Ldim + 0.01;
               Pmax = 0.5/Ldim +0.01; 
            end
            
       k=0;
       
       while (Pend<Pmax)
            %Pend = Pend+dPend;
            Pend = Pend+onecyl;
            k=k+1;
        end
        %cyladd = [cyladd; k*dPend];
        cyladd = [cyladd; k*onecyl];
        
        y0(1:K) = y(end,1:K);
        y0(1+K:2*K) = y(end,K+1:2*K) ;
        y0(2*K+1:3*K) = y(end,2*K+1:3*K);
        y0(3*K+1) = Pend ;
        y0(3*K+2) = y(end,3*K+2);  
        %tspan = [0 tend];
        tout = linspace(tstart,tend,N);
    end


%tout = linspace(te,T,N) ; 
options = odeset('RelTol',1.0e-4,'AbsTol',1.0e-4);
    
%options = odeset('RelTol',1.0e-03,'AbsTol',1.0e-06);

[t,y] = ode15s(@coupledPdeWithdPdt,tout,y0,options); 
%[t,y] = ode15s(@coupledPdeWithdPdt,tout,y0,options); 

%first = [first; te]; 
%The 't' values are given at the rows 

Alamt   = [Alamt; y(:,1:K)]; % This is A from X=0 to X=1 (this is, 0<x<lambda)

phi =  [phi; y(:,2*K+1:3*K)];

P   = [ P; y(:,3*K+1)]; 
 
lam = [lam; y(:,3*K+2)]; 

A = y(:,1:K)./y(:,3*K+2);
th  = [th; y(:,K+1:2*K)./A];

tvec = [tvec; t];
tstart = tend; 
end
A = Alamt./lam; 
%avg = sum(cyladd)/length(first) 
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

uint  = [ u , ...
					uf.*ones(N,K+1) ];
temp = [th, phi];		% complete temperature profile (theta and phi)


xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

%indx = find(t==dvalues(27),1);
%indx = indx+1; 
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
uftval = 0; 
PlottingContoursAlternative
