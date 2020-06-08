% Clear previous files
function SS = ComputeAndSaveSS(K)
%% COMPUTE STEADY STATE
% Define parameters
% Parameters shared with other routines (alternatively you can compute them
% separately with the dimensional parameters) 

global N T L 


        
  % Do you want the Heaviside? (Yes, you do). 
        H=1;

        % Plots for steady state - 1 , no plots for steady state - 0
        plt = 0;
        eps = 1e-4;

        [~, A0steady, ~, th0steadyfull, xsteady] = InitialConditionsSteady(eps,H,plt);
        % A0 and th0 are actually whatever size Matlab needs them to be, as it uses
        % an adaptive mesh. It is our job to interpolate this accordingly. 
        %close all 

       % Find lam0, then resize both sides with an interpolation (only necessary
        % to do this for the steady state, the conditions on the rest are much
        % nicer because of our rescalings 
        
        %xavg   = (xsteady(1:end-1) + xsteady(2:end))/2;
        I =    find(A0steady>0.97, 1,'first');

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

                     
        


%% START FROM HERE WHEN YOU HAVE ALREADY CALCULATED STEADY STATE 

y0(1:K) = A0steady;
y0(1+K:2*K) = A0steady.*th0steady;
y0(2*K+1:3*K) = phi0steady;		
 
y0(3*K+1) = lam0steady; 

% Independent variable for ODE integration 
tout = linspace(0,T,N);

%% ODE integration 
options = odeset('RelTol',1.0e-06,'AbsTol',1.0e-06);

tic
[t,y] = ode15s(@coupledPde,tout,y0); 
toc

A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

phi   = y(:,2*K+1:3*K);

lam = y(:,3*K+1); 

% Save the solution from here and then import into steady state to see if
% it actually converges to a steady state
SS = [A(end,:)'; th(end,:)'; phi(end,:)'; lam(end)]; 
csvwrite('SteadyStateK600.csv', SS); 
end