close all 
clear all

global N K T D uf P0 uft P0t
% We define all of the parameters in an external routine for clarity 
ParametersDefinition

%P0 = c1*D*St
% set to 1 if we want to compare with steady state
st = 1; 

if st == 1 
    % Change filename to match what we want to import 
    data = csvread('SSK300uft1.csv');
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
        %P0   = data(10*K+4);
end

%% START FROM HERE WHEN YOU HAVE ALREADY CALCULATED STEADY STATE 

% Initial conditions 
% 
% A0 = A0steady; 
% th0 = th0steady; 
% th0bot = phi0steady; 

%Q = Bi; % can my code cope with that? Even if the steady state code can't. 
%[Probably not but sort of worth a chance]
y0(1:K) = A0lamt;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = phi0;		
y0(3*K+1) = lam0; 

val = linspace(0.5,1,11)';
val = [val(1:end-1); linspace(1, 5.0, 20)'];

%val = [0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5];
%val = 0.8;
% In this case (unfortunately) I want to control the length of tout

data = [];
omegatmat = []; 
options = odeset('RelTol', 1.0e-4, 'AbsTol', 1.0e-4);

Deltau = 0.5; 
%omega = 2 ;
P0 = 1;
P0t = @(t) P0; 
Tmax = 5; 
%% ODE integration 
for j =1:length(val)
    omega = val(j);
    tout = linspace(0,2*pi*Tmax/omega,N);
    uft = @(t) uf*(1 + Deltau.*sin(omega*t)); % base case 
    % Independent variable for ODE integration 
    %tspan = [0 T];
    tic
    
   % [t,y] = ode15s(@coupledPde,tspan,y0); 
    [t,y] = ode15s(@coupledPdeuft,tout,y0); 
    toc
    %N = length(t);
    Alamt  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)
    lam = y(:,3*K+1); 
    A = Alamt./lam;
    th = y(:,K+1:2*K)./A;
    phi   = y(:,2*K+1:3*K);

    % Saving data
    
    indx = find(tout>=(2*pi*(Tmax-1))/omega,1);
    data = [data, lam(indx:end)];
    omegatmat = [omegatmat, omega*tout(indx:end)']; 
    
    y0(1:K) = Alamt(indx,:);
    y0(1+K:2*K) = A(indx,:).*th(indx,:);
    y0(2*K+1:3*K) = phi(indx,:);		
    y0(3*K+1) = lam(indx); 
end

[r,~] = size(omegatmat); 
valmat = ones(r,1)*val';

min(min(data))
max(max(data))

figure;
contourf(omegatmat-2*pi*(Tmax-1),valmat,data,30,'LineColor','none')
colormap(jet);
xlim([0 2*pi])
% 
 csvwrite('4piminusomegatmatut.csv',4*pi-omegatmat);
 csvwrite('valmatut.csv', valmat);
 csvwrite('lamdataut.csv', data); 
%axis off
return 

print(gcf, '-dpng', '-r600', '-painters', 'ufomegaSweepjet.png')

