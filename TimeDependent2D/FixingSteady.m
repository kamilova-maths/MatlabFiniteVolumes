% Fixing steady state. DELETE after getting it done.

% Clear previous files
close all 
clear all
clc


%% COMPUTE STEADY STATE
% Define parameters
% Parameters shared with other routines (alternatively you can compute them
% separately with the dimensional parameters) 

global Pe Bi tha N K gamma P0 St T L D uf x1 x2 Q P0t

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
T_a=343; 
h = 7; 

%Defining non-dimensional parameters
% Peclet number
Pe = (rho*c*uc*Ldim)/(k);
epsilon=R1/Ldim;
St=(rho*g*Ldim^2)/(uc*mu0);

P0 = (10000*Ldim)/((R1^2)*uc*mu0);
Bi= ((Ldim^2)*h)/(k*R1); 
DeltaT = (Qc*Ldim)/(rho*c*uc);
%DeltaT = (Qc*Ldim)/(Bi*rho*c*uc); 

tha = 0.005; 
D = (R0^2)/(R1^2); 

gamma = 20; 

%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Q = 1;

uf = 1; 
% Pe = 37.8; St = 8.8; P0 =0.7; Bi = 114.3; tha=0.005; D = 0.25; 
% gamma = 30;  x1 = 5/7; x2 = 6.5/7; Q = 1; uf = 1; 

% Calculating the initial conditions as a solution of the steady state
% problem 
% Discretisation in t
N=800; 
% Discretisation in x
K=800;

% end of the domain
T = 2; L=1.5 ; 

% Calculate steady state 
        % Do you want the Heaviside? (Yes, you do). 
H=1;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;
eps = 1e-4;

[Psteady, A0steady, J0steady, th0steadyfull, xsteady] = InitialConditionsSteady(eps,H,plt);
% A0 and th0 are actually whatever size Matlab needs them to be, as it uses
% an adaptive mesh. It is our job to interpolate this accordingly. 
%close all 

% Find lam0, then resize both sides with an interpolation (only necessary
% to do this for the steady state, the conditions on the rest are much
% nicer because of our rescalings 

%xavg   = (xsteady(1:end-1) + xsteady(2:end))/2;
I =    find(A0steady>0.97, 1,'first');
I = I+12; 


lam0steady = xsteady(I);

A0intfull = interp1(xsteady,A0steady,linspace(0,L,K+1)','pchip'); 
A0celfull = (A0intfull(1:end-1) + A0intfull(2:end))/2; 
A0int   = interp1(xsteady(1:I),A0steady(1:I),linspace(0,xsteady(I),K+1)','pchip'); 


th0intfull = interp1(xsteady,th0steadyfull,linspace(0,L,K+1)','pchip'); % pchip and cubic should be exactly the same
th0celfull = (th0intfull(1:end-1) + th0intfull(2:end))/2; 

th0int  = interp1(xsteady(1:I),th0steadyfull(1:I),linspace(0,xsteady(I),K+1)','pchip'); % pchip and cubic should be exactly the same

%         A0int   = interp1(xsteady,A0steady,linspace(0,L,K+1)','pchip'); 
%         th0int  = interp1(xsteady,th0steadyfull,linspace(0,L,K+1)','pchip'); % pchip and cubic should be exactly the same



phi0int = interp1(xsteady(I+1:end),th0steadyfull(I+1:end),linspace(xsteady(I+1),L,K+1)','pchip');
%         
A0cel   = (A0int(1:end-1)+A0int(2:end))/2; % this is  size K x 1  - cells   
th0cel  = (th0int(1:end-1) + th0int(2:end))/2;  % this size K x 1  -cells 
phi0cel = (phi0int(1:end-1) + phi0int(2:end))/2; % this is size K x1 - cells

% Actually, I should just use the full A0 and th0, as long as I
% extract lambda correctly... right? This is all just extra effort

% Define the steady states
A0steady      = A0cel;
th0steady     = th0cel;
phi0steady    = phi0cel; 
u0steady      = 1./A0int;

A0steadyfull = A0celfull; 
u0steadyfull = 1./A0intfull; 
th0steadyfull = th0celfull; 
             
% Compare with new computed state for u
 
ucomp = usolution(A0steady,th0steady,1.0564,1,P0);  
ucomp = [ucomp; uf]; 

plot(linspace(0,1,K+1),u0steady,linspace(0,1,K+1),ucomp)