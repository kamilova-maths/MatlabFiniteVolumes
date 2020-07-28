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
    data = csvread('SSK300.csv');
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

val1 = linspace(0.1,1,2)';
val2 = [val1'; val1']; 
val = val2(:); 

%val = [0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5];
%val = 0.8;
% In this case (unfortunately) I want to control the length of tout
tout = linspace(0,T,N);
data = [];
valvec = []; 
omega = 2 ;
P0 = 1;
P0t = @(t) P0; 
%% ODE integration 
for j =1:length(val)
    Deltau = val(j);
    %options = odeset('RelTol',[1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3],'AbsTol',[1.0e-04, 1.0e-06, 1.0e-06, 1.0e-04]);
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

    % We calculate u with the solution for A, th and lam
    u = zeros(size(A));
    for i=1:N
        u(i,:) = usolutionuft(A(i,:)',th(i,:)',lam(i),1,P0t(t(i)),uft(t(i)));   
    end
    
    if mod(j,2) == 0
        valvec = [valvec, val];
        data = [data, lam];
    end;
    
    indx = N-1;
    y0(1:K) = Alamt(indx,:);
    y0(1+K:2*K) = A(indx,:).*th(indx,:);
    y0(2*K+1:3*K) = phi(indx,:);		
    y0(3*K+1) = lam(indx); 
%     
%     y0(1:K) = A0lamt;
%     y0(1+K:2*K) = A0.*th0;
%     y0(2*K+1:3*K) = phi0;		
%     y0(3*K+1) = lam0; 


end

valmat = ones(N,1)*val1';
omegatmat = tout'.*omega*ones(size(val1')); 
min(min(data(1:indx,:)))
max(max(data(1:indx,:)))

figure;
contourf(tout',val1',data',20,'LineColor','none')
colormap(jet)

figure; 
contourf(omegatmat,valmat,data,20,'LineColor','none')
colormap(jet)
xlim([0 2*pi])
return
indx = find(tout>1,1);
contourf(tout(1:indx)',val1',data(1:indx,:)',20,'LineColor', 'none')
colormap(jet)
min(min(data(1:indx,:)))
max(max(data(1:indx,:)))
csvwrite('P0DeltaParameterSweep.csv',[tout(1:indx)', data(1:indx,:)]); 

axis off
print(gcf, '-dpng', '-r600', '-painters', 'P0DeltaSweep.png')
%  csvwrite('SSData.csv', SS); 
   % disp('Remember to change the name of the file at the end. Include Gamma and K')

% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms
% set to 1 if we want to save data in csv file 
dat = 0; 
sav = 0; % indicator for saving data
P0tval = 1;
%PlottingTimesteps
