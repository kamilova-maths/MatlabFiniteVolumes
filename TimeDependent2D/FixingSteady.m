% Fixing steady state. DELETE after getting it done.

% Clear previous files
close all 
clear all
clc


%% COMPUTE STEADY STATE
% Define parameters
% Parameters shared with other routines (alternatively you can compute them
% separately with the dimensional parameters) 

ParametersDefinition
global K P0 L uf 

% Calculate steady state 
        % Do you want the Heaviside? (Yes, you do). 
H=1;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;
eps = 1e-3;

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