% DATE:     2020 
% DESCR:    TimeDependentMOLIsothermalGridconvergence
%           Main code for time dependent problem in the isothermal case (Gamma = 0).
%           Uses ParametersDefinition, coupledPdeIsothermal, and usolution to solve pdes. 
%           First we calculate an analytical solution, where we choose the
%           Stokes number such that lambda  = 1. Then we calculate u for
%           this Stokes number and compare and print the difference. The
%           discrepancy is attributed to discretisation errors, and there
%           is no way to circumvent it. Then we use either simple or steady
%           conditions to start the code. We impose P0t as a function, but
%           can be chosen to be a constant value if required. The results
%           are plotted at discrete timesteps to show convergence towards
%           steady state. 
% INPUT: 
%           No input variables
%          
% OUTPUT:   Main outcomes: 
%           Acel: The N x 2*K matrix storing the cell values for A 
%           uinterf: The N x (2*K +1) matrix storing all interface values for
%           u
%           
% ADDITIONAL COMMENTS: 
%
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           coupledPdeIsothermal: This is where the isothermal pdes are, that will be then
%           solved with method of lines with ode15s

clear all
close all 

%% FIRST PART: PARAMETER DEFINITION, CALCULATING ANALYTICAL SOLUTION FOR COMPARISON

% Parameters shared with other routines 
global uin P0t

% Define parameters
ParametersDefinition

% We add the heaviside with H=1, and we remove it with H=0. 
H     = 0;
Gamma = 0; 
uin   = @(t) 1/D; % this is what it should be 
%uin = @(t) 3.9946; % this is what it converges to 
%uin = @(t) (3.9946).*(1+sin(2*pi*t)); 

% Plots for steady state - 1 , no plots for steady state - 0
plt   = 0;
P0t =@(t) P0;

%% CALCULATING ANALYTIC SOLUTION FOR CHOSEN PARAMETERS.
% NOTE THAT THE EXACT SOLUTION IS ONLY DEFINED UP TO LAM0. It is assumed
% that A is 1 afterwards and lambda stays constant. 
% Analytic solution for A and P0 for the steady state, with
% constant mu (no temperature) 
lamex = @(St) atan(sqrt(P0^2 + 6*St - 6*D*St)/sqrt(6*St*D-P0^2)).*6./(sqrt(6*St*D-P0^2)) - ...
        6*atan(P0./sqrt(6*St*D-P0^2))./(sqrt(6*St*D-P0^2)); 
prompt = ' Set St to 1  or lam to 1 ? ';
char = input(prompt);
switch char 
    case 'lam'
        newSt = fsolve(@(St) lamex(St)-1, 21);
        St    = newSt;
    case 'St'
        St=1;
end

lam0 = lamex(St);
fac   = 6*St*D-P0^2;
P0bar =  6*atan(P0./sqrt(fac))./(sqrt(fac)); 
% 


%% ODE INTEGRATION
reltol  = 1.0e-6; abstol = 1.0e-8; %% CHANGING TOLERANCE DIDN'T DO ANYTHING
options = odeset('RelTol',reltol,'AbsTol',abstol);


% We can either use simple conditions (which satisfy the boundary
% conditions) to evolve our solution to a steady state, or from a nearby
% steady state, to refine a current understanding of steady state. 

% We can either use tout or tspan, depending on user input
Gridvalues = floor(logspace(1,3,10));
err_it = zeros(length(Gridvalues),3); % each row represents a grid value, each column one of lambda, A and u
err_ex = zeros(length(Gridvalues),3); % error with the exact solution
lam_prev = 0; 
A_prev = zeros(1,K);
u_prev = zeros(1,K);

for j=1:length(Gridvalues)
    K = Gridvalues(j);   
    N = 9*Gridvalues(j);
    
    x     = linspace(0,lam0,K+1);
    Aex   = (fac/(6*St)).*tan(sqrt(fac).*(x+P0bar)./6).^2 - P0.^2./(6*St)+D ; 


% Note that Aex is A evaluated at the nodes. I want to evaluate it at the
% cells, so I have to perform an averaging, of the form 
    
    uf   = 1;
    %Aex  = [Aex'; 1];
    Aex   = (Aex(1:end-1) + Aex(2:end))./2; 
    uex   = 1./Aex; 
    
     %% AWAY FROM STEADY STATE
    lamhat = 0.7;
    Ainitial = (1- D).*linspace(0,lamhat,K)'./lamhat+D; 
    %Ainitial = (Ainitial(1:end-1)+Ainitial(2:end))/2;
    y0(1:K)  = Ainitial*lamhat; 
    y0(1+K)  = lamhat;



    tout  = linspace(0,T,N);
    [t,y] = ode15s(@coupledPdeIsothermal,tout,y0,options);



    lam   = y(:,K+1);
    A     = y(:,1:K)./lam;
    u     = zeros(size(A));

    for i=1:N
        u(i,:) = usolution(A(i,:)',zeros(size(A(i,:)))',lam(i),1,P0t(t(i))); 
    end

    %  We construct u from A 
   % Acel  = [ A, ones(N,K)];		
   % Aint  = ([2*D*ones(N,1) - A(:,1), A ] + [ A, ones(N,1)] )/2;  % A (interfacess)
    % u interfaces
    % uinterf  = [ u , ...
      %      uf.*ones(N,K+1) ];   
    
    err_it(j,1) = max(abs(lam_prev-lam(end)));
    err_it(j,2) = max(abs(max(A_prev)-max(A(end,:))));
    err_it(j,3) = max(abs(max(u_prev)-max(u(end,:)))); 
    
    err_ex(j,1) = max(abs(lam0-lam(end)));
    err_ex(j,2) = max(abs(Aex-A(end,:)));
    err_ex(j,3) = max(abs(uex-u(end,:)));
    lam_prev = lam(end);
    A_prev = A(end,:);  
    u_prev = u(end,:);
end
csvwrite('TextFiles/err_ex.csv',[(1./(Gridvalues.*9.*Gridvalues))',err_ex]);
return    
% rewrite below into external routine PlottingTimesteps
%% PLOTTING THINGS - Feel free to comment it out as necessary

   
% We define all the x vectors that we might use.

% X values evaluated from 0 to lambda at the interfaces
xint       = linspace(0,1,K+1)';
% X values evaluated from 0 to lambda at the cells
xcel       = linspace(xint(2)/2,1-xint(2)/2,K)';
% x values from both sides, evaluated at the cells
xvector1   = [xcel*lam(end);lam(end) + xcel*(L-lam(1))];

% steady state x values, evaluated at interfaces
xsteadyint = linspace(0,lam0,K)';

% steady state x values, evaluated at the cells
xsteadycel = linspace(xsteadyint(2)/2,1-xsteadyint(2)/2,K)'; 

% Number of timesteps to plot
numel = 10; 

% Plotting A and comparing wiht steady state 
figure; 

% Indices extraction so that plotting data in pgfplots does not take an
% enormous amount of time, and does not produce incredibly slow and large
% figures

Kindices = [1, 2:1:2*K]; 

Adata    = [xvector1(Kindices), Acel(1,Kindices)'];
% We plot the initial condition in dashed lines

plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], Acel(1,:)', '--');
hold on

% We plot the numel timesteps for A, and save the data 
for i = floor(N/numel):floor((N/numel)):N
    xvector = [xcel*lam(i);lam(i) + xcel*(L-lam(i))];
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], Acel(i,:)');
    Adata   = [Adata, [xvector(Kindices), Acel(i,Kindices)']];
end
%
hold on 

% We plot the steady state 
plot(xsteadycel,A0,'--')

% Plotting u 

figure; 
% We plot the initial condition 

plot([xint*lam(1);lam(1) + xint(2:end)*(L-lam(1))], uinterf(1,:)', '--');
udata = [xvector1(Kindices), uinterf(1,Kindices)'];
hold on
for i = floor(N/numel):floor((N/numel)):N
    xvector = [xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))];
	plot([xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))], uinterf(i,:)');
    udata   = [udata, [xvector(Kindices), uinterf(i,Kindices)']];
end
% We plot the steady state 
plot(linspace(0,lam0,K),u0,'--')

figure;
plot(t,lam)
hold on 
plot(t,P0t(t))
xlabel('$t$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

%
prompt = 'Do you want to save the data? (yes == 1 ) ';

sav  = input(prompt); % indicator for saving .csv files. sav == 1 saves files. 

if sav == 1
    csvwrite('TextFiles/Awithlambda1constantmu.csv',Adata); 
    csvwrite('TextFiles/ASteadywithlambda1constantmu.csv', [xsteadycel([1, 5:5:K]), A0([1, 5:5:K])])
    csvwrite('TextFiles/uwithlambda1constantmu.csv',udata); 
    csvwrite('TextFiles/uSteadywithlambda1constantmu.csv', [xsteadyint([1, 5:5:K]), u0([1, 5:5:K])])
    csvwrite('TextFiles/lambdaconstantmup4.csv', [t(1:20:end),lam(1:20:end)])
end
