% Data to use 

rho= 1.8*10^3; %Bergstrom ; 
g = 10; 
c= 900; % Fitt and Howell
L=7; % Temperature Profiles in Soderberg Electrodes
uc = 10^-5; %Bergstrom approximation
R0=0.5;
R1=1; 
k = 3; 
Qc=15000; %Taken very vaguely from Temperature profiles in Soderberg electrodes. 
mu0=10^10; % given by Bjornar at a reference temperature
T_a=343; 
h = 7; 

%Defining non-dimensional parameters
% Peclet number
Pe = (rho*c*uc*L)/(k);
epsilon=R1/L;
St=(rho*g*L^2)/(uc*mu0);
P0 = (10000*L)/((R1^2)*uc*mu0);

DeltaT = (Qc*L)/(rho*c*uc);
Bi= ((L^2)*h)/(k*R1); 
tha = 0.005; 
D = (R0^2)/(R1^2); 


%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
eps = 1e-2;

% In order to set constant viscosity, we set gamma=0. 
gamma = 0; 

% We compute the analytical solution. 
% We solve for the integration constants, with initial guess given by
% Mathematica code, in order to (correctly) satisfy the boundary
% conditions.

% This depends highly on the initial guesses here, it is important to
% choose the ones guided by the plot. (Maybe justify with MATLAB plot?)

% THESE ARE THE CHOICES OF CONSTANTS THAT MATCH THE STEADY STATE SOLUTION I
% HAVE. 
c0 = [-6.5,-4]; 


c = fsolve(@(c) AnalyticalSolutionConstantViscosity(c,St,P0,D),c0); 

% Then the analytical solution for P and A is given by these formulas

eta= 1/3; 

Pan = @(x) (1/sqrt(eta)).*( sqrt(2).* sqrt (eta .*c(1) + eta * c(1)* ...
    (-1 + tanh( 0.5 * ( -sqrt(2).*x.*eta.*sqrt(c(1)) - sqrt(eta) .*sqrt(c(1)) .*c(2))).^2)));

Aan = @(x) (eta.*c(1)./St).*(-1 +tanh(0.5.* (-sqrt(2).*x.*eta.*sqrt(c(1))-sqrt(eta).*sqrt(c(1)).*c(2) )).^2 );

% Now we attempt to plot what is happening with the constants and why we
% have so many different intersections. Is it possible to find a
% combination that will give us the steady state produced by the
% time-dependent code? Intuition says no. 

F1 = @(c1,c2) sqrt(2).* sqrt (c1*(tanh(0.5.*sqrt((1/3).*c1).*c2)).^2)-P0; 
F2 = @(c1,c2) St.*D+(1/3).*c1.*(sech(0.5.*(sqrt((1/3).*c1).*c2))).^2;

fimplicit(F1,[-10 0 -10 10])
hold on 
fimplicit(F2,[-10 0 -10 10])

% Now we compute the Steady state with our original steady state code, for
% constant viscosity

% Number of x nodes
N=100 ; 
L= 1.4529; 
[Pnum, Anum, Jnum, Thnum] = InitialConditionsSteady(N,gamma,Q,x1,x2,eps,St,tha,Bi,Pe,P0,R0,L);

%x = linspace(0,L,N); 

% figure(1)
% plot(x,Anum)
% hold on 
% plot(x,Aan(x),'--')
% 
% 
% figure(2)
% 
% plot(x,Pnum)
% hold on 
% plot(x,Pan(x),'--')



% PERFECT MATCH - This does NOT match the solution I get from the time
% dependent problem


% We definitely use the correct steady state to solve for the
% time-dependent problem. Now we compare with the code for the
% time-dependent problem

K=2500; 
% end of the domain
T = 1; 
th0=Thnum';
A0 = Anum';


[ th, A, u, x, t ] = TimeDependentFDfull_v3( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N );

