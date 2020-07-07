% Clear previous files
clear all 
close all 
%clc

%% MAIN CODE FOR 'EASY' VERSION OF EXPERIMENTS

% These are a series of tests, which currently ALL fail. I am absolutely
% all out of ideas for sure. 
% First we define the parameters (feel free to ignore it) 
% Then we solve the problem analytically (feel free to ignore that too) 
% NOTE: THE ANALYTICAL SOLUTION IS CHOSEN SUCH THAT LAMBDA = 1 

% TEST 1: Early comparison, we compare u obtained analytically with the u
% obtained with my code, using the A obtained analytically as input and
% lambda = 1

% TEST 2: initial = 'steady', Pvalue = 'constant' 
% We use the steady state (analytical) solution as initial condition, and P
% is just a constant (no dP/dt equation). This should converge to lambda =
% 1 forever, and A and u remain unchanged. In reality, u is changed and
% lambda is not 1, (although it's close to 1). I still consider this a
% failed test because this is the easiest possible condition. 

% TEST 3: initial = 'steady', Pvalue = 'dPdt'
% We use the steady state (analytica) solution as initial condition, and P
% is the solution of a first order ODE. The equation is 

% dP/dt = D*St(uin(t)-u(1))  

% meaning that, if uin = u(1), I should have a steady state. This is why
% it's important that I understand what u(1) is ... 
% This test runs with uin = @(t) 1/D, which should be the right condition.
% Alternatively you can make uin to be whatever it converges to and re-run
% it.THIS ALSO DOES NOT WORK 

% Further tests are done with initial = 'simple', where we use simple
% initial conditions and see if it converges to the steady state. 


%



%% FIRST PART: PARAMETER DEFINITION, CALCULATING ANALYTICAL SOLUTION FOR COMPARISON

% Parameters shared with other routines 
global Pe Bi N K Gamma P0 St T L D uf uin

% Define parameters
% Data to use 

rho= 1.8*10^3; %Bergstrom ; 
g = 10; 
c= 900; % Fitt and Howell
Ldim=5; % Temperature Profiles in Soderberg Electrodes
uc = 10^-5; %Bergstrom approximation
R0=0.5;
R1=1; 
k = 3; 
Qc=15000; %Taken very vaguely from Temperature profiles in Soderberg electrodes. 
%Qc = 3.5607e+03;
%Qc = 2000/L; 
%mu0=10^10; % given by Bjornar at a reference temperature
mu0 = 4E10; 
T_in= 20 + 273; 
T_a = 80 + 273; 
h = 7; 

%Defining non-dimensional parameters
% Peclet number
Pe  = (rho*c*uc*Ldim)/(k);
epsilon = R1/Ldim;
%St=(rho*g*Ldim^2)/(uc*mu0);
St = 27; 
%P0 = (10000*Ldim)/((R1^2)*uc*mu0);
%P0=0;
Bi= ((Ldim^2)*h)/(k*R1); 

P0 = 0.1;
%tha = 0.005; 
D = (R0^2)/(R1^2); 

% Calculating the initial conditions as a solution of the steady state
% problem 
N=2000; K=400;
% end of the domain
T = 100; L= 1 ;


% We add the heaviside with H=1, and we remove it with H=0. 
H     = 0;
Gamma = 0; 
%uin  = @(t) 1- 1/D; % this is what it should be 
%uin = @(t) 0.2; 
%uin = @(t)  3.9953; % this is what it converges to 
uin = @(t) (1/D).*(1+0.5*sin(2*pi*t)); 
%uin = @(t)0.5*(1+0.5*sin(2*pi*t));
%uin = @(t) (0.2).*(1+0.5*sin(2*pi*t)); 

% Plots for steady state - 1 , no plots for steady state - 0
plt   = 0;
%% CALCULATING ANALYTIC SOLUTION FOR CHOSEN PARAMETERS.
% Analytic solution for A and P (and u?) 0 for the steady state, with
% constant mu (no temperature) 
lamex = @(St) atan(sqrt(P0^2 + 6*St - 6*D*St)/sqrt(6*St*D-P0^2)).*6./(sqrt(6*St*D-P0^2)) - ...
    6*atan(P0./sqrt(6*St*D-P0^2))./(sqrt(6*St*D-P0^2)); 

newSt = fsolve(@(St) lamex(St)-1, 21);
St    = newSt;

x     = linspace(0,L,K);
fac   = 6*St*D-P0^2;
P0bar =  6*atan(P0./sqrt(fac))./(sqrt(fac)); 
Aex   = (fac/(6*St)).*tan(sqrt(fac).*(x+P0bar)./6).^2 - P0.^2./(6*St)+D ; 



lam0 = lamex(St); 
% Note that Aex is A evaluated at the nodes. I want to evaluate it at the
% cells, so I have to perform an averaging, of the form 
u0 = 1./Aex'; 
uf = 1;
Aex = [Aex'; 1];
%Aex = [D; Aex'];
A0 = (Aex(1:end-1) + Aex(2:end))./2; 


%% FIRST, EARLY COMPARISON: I am not sure this is worth doing 
% I'm not sure if the problem is here, or if this matters. 
% Here I compare the calculated u with the exact u, inputting the
% analytical value of A
ucalc = usolution(A0,zeros(size(A0)),lam0,L,P0);
max(abs(ucalc-u0))

%return

%% USING TIME DEPENDENT SOLVER 
% Initial conditions that satisfy the boundary conditions but are not the
% analytical solution ( apart from lambda)

% I don't solve for theta as this is the constant viscosity case, but I
% still input it as a zero vector in order to keep the structure of the
% numerical method
th0 = zeros(size(A0)); 

% dirichlet condition for u at x=1

initial = 'steady';
switch initial
    case 'simple'
        %% AWAY FROM STEADY STATE
        Ainitial = (1- D).*linspace(0,1,K+1)'+D; 
        Ainitial = (Ainitial(1:end-1)+Ainitial(2:end))/2; 
        % We now check the return to steady state 
        y0(1:K) = Ainitial; 

        y0(1+K) = 0.8;

    case 'steady'
        %% STEADY STATE
        y0(1:K) = A0; 
        y0(1+K) = lam0;

end

%% ODE INTEGRATION

tout = linspace(0,T,N);

reltol = 1.0e-4; abstol = 1.0e-8; %% CHANGING TOLERANCE DIDN'T DO ANYTHING
options = odeset('RelTol',reltol,'AbsTol',abstol);
%Pvalue = 'constant';
Pvalue = 'dPdt';

switch Pvalue
    case 'constant'
        % no need to define y0(2+K)
        [t,y] = ode15s(@coupledPdeNoTemp,tout,y0);
        A  = y(:,1:K);
        lam   = y(:,1+K);
        u = zeros(size(A));
        %umat(1,:) = u0; 
        for i=1:N
            u(i,:) = usolution(A(i,:)',zeros(size(A(i,:)))',lam(i),L,P0); 
        end
        
    case 'dPdt'
        y0(2+K) = P0;
        [t,y] = ode15s(@coupledPdeNoTempdPdt,tout,y0); 
        A  = y(:,1:K);
        lam   = y(:,1+K);
        P = y(:,2+K);
        u = zeros(size(A));
        %umat(1,:) = u0; 
        for i=1:N
            u(i,:) = usolution(A(i,:)',zeros(size(A(i,:)))',lam(i),L,P(i)); 
        end
end

%% PLOTTING THINGS - Feel free to comment it out as necessary 

%  We construct u from A 
Acel = [ A, ones(N,K)];		
Aint = ([2*D*ones(N,1) - A(:,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)

   
% We rescale it as 

uint  = [ u , ...
		uf.*ones(N,K+1) ];   
   
dx = 1/K;

x = (0:dx:L)';

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
xvector1 = [xcel*lam(1);lam(1) + xcel*(L-lam(1))];
figure; 
numel = 10; 
Kindices = [1, 5:5:2*K]; 

%Adata = [xvector1(Kindices), Acel(1,Kindices)'];
plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], Acel(1,:)', '--');
	hold on
for i = N/numel:(N/numel):N
    xvector = [xcel*lam(i);lam(i) + xcel*(L-lam(i))];
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], Acel(i,:)');
 %   Adata = [Adata, [xvector(Kindices), Acel(i,Kindices)']];
end
%csvwrite('Awithlambda1constantmu.csv',Adata); 
hold on 
xsteadyint = linspace(0,lam0,K)';
xsteadycel = linspace(xsteadyint(2)/2,1-xsteadyint(2)/2,K)'; 
plot(xsteadycel,A0,'--')
%csvwrite('ASteadywithlambda1constantmu.csv', [xsteadycel([1, 5:5:K]), A0([1, 5:5:K])])

% Plot it the right way, K
figure; 
plot([xint*lam(1);lam(1) + xint(2:end)*(L-lam(1))], uint(1,:)', '--');
%udata = [xvector1(Kindices), uint(1,Kindices)'];
hold on
for i = N/numel:(N/numel):N
    xvector = [xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))];
	plot([xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))], uint(i,:)');
 %   udata = [udata, [xvector(Kindices), uint(i,Kindices)']];
end
plot(linspace(0,lam0,K),u0,'--')
%csvwrite('uwithlambda1constantmu.csv',udata); 
%csvwrite('uSteadywithlambda1constantmu.csv', [xsteadyint([1, 5:5:K]), u0([1, 5:5:K])])
figure;
plot(t,lam)
if strcmp(Pvalue,'dPdt') == 1
    hold on 
    plot(t,P)
end
xlabel('$t$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

%csvwrite('lambdaconstantmu.csv', [t,lam])

return 