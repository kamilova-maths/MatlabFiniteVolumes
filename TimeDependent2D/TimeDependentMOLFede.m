% Clear previous files
close all 
clc


%% COMPUTE STEADY STATE
% Define parameters
% Parameters shared with other routines 
global Pe Bi tha N K gamma P0 St T L D uf x1 x2 Q
Pe = 37.8; St = 8.8; P0 =0.7; Bi = 114.3; tha=0.005; D = 0.25; 
gamma = 20;  x1 = 5/7; x2 = 6.5/7; Q = 1; uf = 1; 

% Calculating the initial conditions as a solution of the steady state
% problem 
% Discretisation in t
N=800; 
% Discretisation in x
K=300;

% end of the domain
T = 1; L=1.5 ;

 
%Initial conditions with K = 300 
data = csvread('InitialConditionsK300.csv');

A0 = data(1:K); 

th0=data(K+1:2*K); 
th0bot = data(2*K+1:3*K); % We resize to Kx1 and flip the bottom part for theta,
                        % to match the requirements from coupledPde
lam0= data(3*K+1) ; 
Fk0 =  0; % we have to set this to zero, and our equations must be such that this is conserved

%% START FROM HERE WHEN YOU HAVE ALREADY CALCULATED STEADY STATE 

% y0(1:K) = A0;
% y0(1+K:2*K) = A0.*th0;
% y0(2*K+1:3*K) = th0bot;
% y0(3*K+1) = lam0; 
% % If you want the option without mass matrix, just comment the next line,
% % and use corresponding options
% y0(3*K+2) = Fk0;  % This is the flux equation that puts phi and theta together
% % Mass matrix 
% M = speye(3*K+2,3*K+2); 
% M(end,end) = 0;  % This sets the algebraic part of our equation

y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
% th0bot(end) is == th0(end), so drop that extra variable. Also, flip it
y0(2*K+1:3*K-1) = flip(th0bot(1:end-1));
y0(3*K) = lam0; 


% Independent variable for ODE integration 
tout = linspace(0,T,N);

%% ODE integration 
% OPTIONS WITHOUT MASS MATRIX 

options = odeset('RelTol',1.0e-04,'AbsTol',1.0e-04);

% OPTIONS WITH MASS MATRIX

% optionsthtmp = odeset('Mass',M,'RelTol',1e-4,'AbsTol',[1e-4 1e-6 1e-6 1e-4 1e-6]);

tic
[t,y] = ode15s(@coupledPdeFede,tout,y0); 
toc

A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

phi   = y(:,2*K+1:3*K-1);

lam = y(:,3*K); 

% FK = y(:,3*K+2); %
% We calculate u with the solutuon for A, th and lam
u = zeros(size(A));
for i=1:N
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i),1);   
end

% We add the Dirichlet boundary conditions 
thtop = [ zeros(N,1), ...
       th         ];
  
Atop  = [ D*ones(N,1), ...
       A           ];
   
tmpA =  (Atop + [Atop(:,2:end), ones(N,1)])/2; % extract A at the edges 
Abot = ones(N,K);  % At the bottom, A is just one

Afull =    [tmpA, Abot];  
ufull  = [ u , ...
      uf.*ones(N,1) ];
   
dx = 1/K;


% We rescale X and Xbar in order to plot. Note that at the top, where we
% use X, we have K+1 terms, whereas at the bottom, where we use Xbar, we
% have K terms

x = (0:dx:1)';
[X, T1] = meshgrid(x,t);
Xresc1 = lam.*X; 
xbar = linspace(0,1,K-1)'; 
[X, T2] = meshgrid(xbar(2:end),t);
Xresc2 = X.*(L-lam)+lam;

% PLOTTING THETA
% We have to flip phi because Xbar goes from 1 to 0.

figure; 
plot(Xresc1(1,:),thtop(1,:))
hold on 
plot(Xresc2(1,:),phi(1,2:end))
hold on 
numel = 10; 
for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),thtop(i,:))
    hold on
    plot(Xresc2(i,:),phi(i,2:end)) 
end
set(gca,'TickLabelInterpreter','latex','fontsize',13)

return 
% IF you want to see the other solutions, remove the return. 
% PLOTTING A
figure; 
plot(Xresc1(1,:),tmpA(1,:))
hold on 
plot(Xresc2(1,:),Abot(1,:))


for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),tmpA(i,:))
    hold on
    plot(Xresc2(i,:),Abot(i,:))
end



% PLOTTING U
figure; 
plot(Xresc1(1,:),ufull(1,:))
hold on 
plot(Xresc2(1,:),Abot(1,:))


for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),ufull(i,:))
    hold on
    plot(Xresc2(i,:),Abot(i,:))
    
end

