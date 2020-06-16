global Pe Bi tha N K Gamma P0 St T L D uf x1 x2 Q P0t

rho= 1.8*10^3; %Bergstrom ; 
g = 10; 
c= 900; % Fitt and Howell
Ldim=7; % Temperature Profiles in Soderberg Electrodes
uc = 10^-5; %Bergstrom approximation
R0=0.5;
R1=1; 
k = 3; 
Qc=15000; %Taken very vaguely from Temperature profiles in Soderberg electrodes. 
mu0=10^10; % given by Bjornar at a reference temperature
T_in= 20 + 273; 
T_a = 80 + 273; 
h = 7; 

%Defining non-dimensional parameters
% Peclet number
Pe = (rho*c*uc*Ldim)/(k);
epsilon=R1/Ldim;
St=(rho*g*Ldim^2)/(uc*mu0);
%St = 10; 
P0 = (10000*Ldim)/((R1^2)*uc*mu0);
Bi= ((Ldim^2)*h)/(k*R1); 


%tha = 0.005; 
D = (R0^2)/(R1^2); 


%This is the area of the clamps, taken from Temperature profiles ... 
x1dim = 5;
x2dim = 6.5; 
x1 = x1dim/Ldim;
x2 = x2dim/Ldim;
Q0 = Qc*(x2dim - x1dim)/Ldim ; 
eps = 1e-4;
%DeltaT = (Qc*Ldim)/(rho*c*uc);
DeltaT =  (Q0*Ldim)/(rho*c*uc); 
gammaBar = 0.069; 
Gamma = gammaBar*DeltaT; 
Q = 1/(x2-x1);

tha = (T_a- T_in)/DeltaT ; 
uf = 1; 
% Pe = 37.8; St = 8.8; P0 =0.7; Bi = 114.3; tha=0.005; D = 0.25; 
% gamma = 30;  x1 = 5/7; x2 = 6.5/7; Q = 1; uf = 1; 

% Calculating the initial conditions as a solution of the steady state
% problem 
% Discretisation in t
N=800; 
% Discretisation in x
K=300;

P0t = @(t)P0; 
% end of the domain
T = 5; L=1.5 ; 