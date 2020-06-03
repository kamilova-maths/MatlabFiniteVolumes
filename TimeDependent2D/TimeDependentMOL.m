% Clear previous files
close all 
clear all
clc


%% COMPUTE STEADY STATE
% Define parameters
% Parameters shared with other routines (alternatively you can compute them
% separately with the dimensional parameters) 

global Pe Bi tha N K gamma P0 St T L D uf x1 x2 Q

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
eps = 1e-4;

uf = 1; 
% Pe = 37.8; St = 8.8; P0 =0.7; Bi = 114.3; tha=0.005; D = 0.25; 
% gamma = 30;  x1 = 5/7; x2 = 6.5/7; Q = 1; uf = 1; 

% Calculating the initial conditions as a solution of the steady state
% problem 
% Discretisation in t
N=800; 
% Discretisation in x
K=300;

% end of the domain
T = 5; L=1.5 ; 

method = 'calculate';
switch method
% Import steady state
    case 'import'
        data = csvread('InitialConditionsK300.csv');

        A0steady = data(1:K); 

        th0steady=data(K+1:2*K); 
        th0botsteady = flip(data(2*K+1:3*K));
        lam0steady = data(3*K+1) ; 
        u0steady = 1./A0steady; % This is not quite accurate, take with a grain of salt. 
        % It's not used for any calculations but it is plotted as a
        % reference when we're too lazy to calculate it again . 
    case 'calculate'
        
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

        lam0steady = xsteady(I);

%         A0int   = interp1(xsteady(1:I),A0steady(1:I),linspace(0,xsteady(I),K+1)','pchip'); 
%         th0int  = interp1(xsteady(1:I),th0steadyfull(1:I),linspace(0,xsteady(I),K+1)','pchip'); % pchip and cubic should be exactly the same

        A0int   = interp1(xsteady,A0steady,linspace(0,L,K+1)','pchip'); 
        th0int  = interp1(xsteady,th0steadyfull,linspace(0,L,K+1)','pchip'); % pchip and cubic should be exactly the same

        u0steady = 1./A0int; 
        
%         phi0int = interp1(xsteady(I+1:end),th0steadyfull(I+1:end),linspace(xsteady(I+1),L,K+1)','pchip');
%         
        A0cel   = (A0int(1:end-1)+A0int(2:end))/2; % this is  size K x 1  - cells   
        th0cel  = (th0int(1:end-1) + th0int(2:end))/2;  % this size K x 1  -cells 
        % phi0cel = (phi0int(1:end-1) + phi0int(2:end))/2; % this is size K x1 - cells
        
        % Actually, I should just use the full A0 and th0, as long as I
        % extract lambda correctly... right? This is all just extra effort
        
        % Define the steady states
        A0steady      = A0cel;
        th0steady     = th0cel;

        
        
end

%% START FROM HERE WHEN YOU HAVE ALREADY CALCULATED STEADY STATE 

% Initial conditions 

A0 = (1- 1e-5 -D).*linspace(0,1,K)'+D; 
th0 = zeros(K,1); 
th0bot = zeros(K,1); 

%Q = Bi; % can my code cope with that? Even if the steady state code can't. 
%[Probably not but sort of worth a chance]
y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = flip(th0bot);		%must flip since th0bot is stored in reverse order
lil = 0.02; 
y0(3*K+1) = 1.0564+lil; 


% Independent variable for ODE integration 
tout = linspace(0,T,N);

%% ODE integration 
options = odeset('RelTol',1.0e-04,'AbsTol',1.0e-04);

tic
[t,y] = ode15s(@coupledPde,tout,y0); 
toc

A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

phi   = y(:,2*K+1:3*K);

lam = y(:,3*K+1); 

% We calculate u with the solution for A, th and lam
u = zeros(size(A));
for i=1:N
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i),1);   
end

% We add the Dirichlet boundary conditions 
thtop = [ zeros(N,1), ...
       th         ];

Acel = [ A, ones(N,K)];														% A (cell values)
Aint = ([D*ones(N,1), A ] + [ A, ones(N,1)] )/2;  % A (interfaces)

uint  = [ u , ...
					uf.*ones(N,K+1) ];
temp = [th, phi];		% complete temperature profile (theta and phi)
				
% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

% PLOTTING THETA

figure;
numel = 10; 
Kindices = [1, 5:5:2*K]; 
xvector1 = [xcel*lam(1);lam(1) + xcel*(L-lam(1))];
plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], temp(1,:)', '--');
thetadata = [xvector1(Kindices), temp(1,Kindices)'];
hold on

for i = N/numel:(N/numel):N
    xvector = [xcel*lam(i);lam(i) + xcel*(L-lam(i))];
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], temp(i,:)');
    thetadata = [thetadata, [xvector(Kindices), temp(i,Kindices)']];
end
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
hold on 
xcelfull= linspace(0,L,K); 
plot(xcelfull,th0steady,'--')
% thsteady = [th0steady; th0botsteady];
csvwrite('ThetaDiscreteTimesteps.csv',thetadata); 

thsteady = [xcelfull([1, 5:5:K])', th0steady([1, 5:5:K])]; 
csvwrite('Thetasteady.csv',thsteady);

title('Temperature')
% Save data to file

% return 
% IF you want to see the other solutions, remove the return. 
% PLOTTING A
figure; 
Adata = [xvector1(Kindices), Acel(1,Kindices)'];
plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], Acel(1,:)', '--');
	hold on
for i = N/numel:(N/numel):N
    xvector = [xcel*lam(i);lam(i) + xcel*(L-lam(i))];
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], Acel(i,:)');
    Adata = [Adata, [xvector(Kindices), Acel(i,Kindices)']];
end
csvwrite('ADiscreteTimesteps.csv',Adata); 

plot(xcelfull,A0steady,'--')
Asteady = [xcelfull([1, 5:5:K])', A0steady([1, 5:5:K])]; 
csvwrite('Asteady.csv',Asteady);

title('Area')
% Save data to file


% PLOTTING U
figure; 
plot([xint*lam(1);lam(1) + xint(2:end)*(L-lam(1))], uint(1,:)', '--');
udata = [xvector1(Kindices), uint(1,Kindices)'];
hold on
for i = N/numel:(N/numel):N
    xvector = [xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))];
	plot([xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))], uint(i,:)');
    udata = [udata, [xvector(Kindices), uint(i,Kindices)']];
end
xintfull = linspace(0,L,K+1); 
plot(xintfull,u0steady,'--')
%usteady = [u0steady(1:end-1); ones(K,1)];
csvwrite('uDiscreteTimesteps.csv',udata); 

usteady = [xintfull', u0steady]; 
csvwrite('usteady.csv', usteady); 

title('Velocity')
% Save data to file


% PLOTTING lambda
figure; 
plot(t, lam);
csvwrite('lam.csv',[[t(1), lam(1)]; [t(5:5:end), lam(5:5:end)]]);
hold on
plot([t(1),t(end)], [lam(1),lam(1)]);
title('lambda')
xlabel('t')

% video version
figure('units','normalized','outerposition',[0 0 0.25 1])
for i = N/numel:(N/numel):N
	subplot(3,1,1)
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], temp(i,:)'), axis([0 L 0 (max(max(temp))+0.1)]);
	hold on
	plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], temp(1,:)', '--'), axis([0 L 0 (max(max(temp))+0.1)]);
	plot([lam(i),lam(i)], [0,max(max(temp))+0.1]);
	hold off
	title(strcat('T=',num2str(i*T/N)))
	ylabel('Temperature')
	subplot(3,1,2)
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], Acel(i,:)'), axis([0 L 0 (max(max(Acel))+0.1)]);
	hold on
	plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], Acel(1,:)', '--'), axis([0 L 0 (max(max(Acel))+0.1)]);
	hold off
	ylabel('Area')
  subplot(3,1,3)
	plot([xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))], uint(i,:)'), axis([0 L 0 (max(max(uint))+0.1)]);
	hold on
	plot([xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))], uint(i,:)', '--'), axis([0 L 0 (max(max(uint))+0.1)]);
	hold off
	ylabel('Velocity')
	pause(T/N)
    %pause(1)
end


