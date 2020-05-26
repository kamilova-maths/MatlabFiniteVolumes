% Clear previous files
clear all 
clc

% Parameters shared with other routines (I'm not very sure about this,but
% let's give it a try)
global Pe Bi tha N K gamma P0 St T L D uf x1 x2 Q

% Define parameters
% Data to use 

rho= 1.8*10^3; %Bergstrom ; 
g = 10; 
c= 900; % Fitt and Howell
Ld=7; % Temperature Profiles in Soderberg Electrodes
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
Pe = (rho*c*uc*Ld)/(k);
epsilon=R1/Ld;
St=(rho*g*Ld^2)/(uc*mu0);
P0 = (10000*Ld)/((R1^2)*uc*mu0);

DeltaT = (Qc*Ld)/(rho*c*uc);
Bi= ((Ld^2)*h)/(k*R1); 
tha = 0.005; 
D = (R0^2)/(R1^2); 

gamma = 20; 

%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
eps = 1e-4;

% Calculating the initial conditions as a solution of the steady state
% problem 
N=1000; K=600;
% end of the domain
T = 1; L=1.5 ;
dx = L/K;
uf = 1; 

% We add the heaviside with H=1, and we remove it with H=0. 
H=1;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;

% We try with tha=0
[P, A0, J0, th0, xsteady] = InitialConditionsSteady(eps,H,plt);
% A0 and th0 are actually whatever size Matlab needs them to be, as it uses
% an adaptive mesh. It is our job to interpolate this accordingly. 
%close all 



% Find lam0, then resize both sides with an interpolation (only necessary
% to do this for the steady state, the conditions on the rest are much
% nicer because of our rescalings 
A0avg  = (A0(1:end-1)+A0(2:end))/2;   
th0avg = (th0(1:end-1) + th0(2:end))/2; 
xavg   = (xsteady(1:end-1) + xsteady(2:end))/2;
I =    find(A0avg>0.99999,1,'first');

lam0 = xavg(I);

A0top = interp1(xavg(1:I),A0avg(1:I),linspace(0,xavg(I),K),'pchip'); 
th0top = interp1(xavg(1:I),th0avg(1:I),linspace(0,xavg(I),K),'pchip'); % pchip and cubic should be exactly the same

%A0bot = ones(size(A0top));
th0bot = interp1(xavg(I:end),th0avg(I:end),linspace(xavg(I),L,K),'pchip');
 
%Initial conditions given by A0, th0
% A0 = A0top'; 
xtemp = linspace(0,1,K);
A0=(A0top + [0, sin(xtemp(1:end-2)), 0])'; %small perturbation in A

%A0  = A1(2:end,end)+[0; sin(x(1:end-2)');0];
A0(end) = 1; 
th0=th0top'; 
th0bot = flip(th0bot)'; % We resize to Kx1 and flip the bottom part for theta,
                        % to match the requirements from coupledPde


y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = lam0.*ones(size(A0));  

y0(3*K+1:4*K) = th0bot;

% Independent variable for ODE integration 
tout = linspace(0,T,N);

% ODE integration 
reltol = 1.0e-04; abstol = 1.0e-04;
options = odeset('RelTol',reltol,'AbsTol',abstol);
%[t,y] = ode15s(@coupledPdeNoTemp,tout,y0);
tic
[t,y] = ode15s(@coupledPde,tout,y0); 
% [t,y] = ode15s(@coupledPde2,tout,y0); 
toc
%pause
A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

lam   = y(:,2*K+1:3*K);

phi = y(:,3*K+1:4*K); 

% This is the solution to u at the top. 
u = zeros(size(A));
for i=1:N
    % Here A(1,:) is A0, so this should give me the u used for the
    % calculations inside coupledPde.m
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i,end),1);   
end
disp('Completed Round 1')

thtop = [ zeros(N,1), ...
       th         ];
  
Atop  = [ D*ones(N,1), ...
       A           ];
   
tmpA =  (Atop + [Atop(:,2:end), ones(N,1)])/2; % extract A at the edges 
Abot = ones(N,K); 
Afull =    [tmpA, Abot];  
ufull  = [ u , ...
      uf.*ones(N,1) ];
   
dx = 1/K;
% we define a vector for lambda(t) [so that I don't get confused with
% matrices]. Recall that lambda only depends on t, so it is constant for
% each x

lamt = lam(:,end); % size  N x 1 [column vector]

x = (0:dx:1)';
[X, T1] = meshgrid(x,t);
Xresc1 = lamt.*X; 
xbar = linspace(1,0,K)'; 
[X, T2] = meshgrid(xbar,t);
Xresc2 = -X.*(L-lamt)+L;

% PLOTTING THETA
numel=10;
datamat1 = [[Xresc1(1,1:5:end)'; Xresc2(1,1:5:end)'], [thtop(1,1:5:end)'; flip(phi(1,1:5:end))']];

figure; 
plot(Xresc1(1,:),thtop(1,:))
hold on 
plot(Xresc2(1,:),flip(phi(1,:)))
hold on 
for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),thtop(i,:))
    hold on
    plot(Xresc2(i,:),flip(phi(i,:))) 
   % pause
    datamat1 = [datamat1, [[Xresc1(i,1:5:end)';Xresc2(i,1:5:end)'], [thtop(i,1:5:end)'; flip(phi(i,1:5:end))']]];
end
set(gca,'TickLabelInterpreter','latex','fontsize',13)
top = max(max(abs(th(1,:)-th(end,:))))/max(max(th))
bottom = max(max(abs(phi(1,:)-phi(end,:))))/max(max(phi))
csvwrite('ThetaDiscreteTimestepsgamma20.csv',datamat1); 

% PLOTTING A
figure; 
plot(Xresc1(1,:),tmpA(1,:))
hold on 
plot(Xresc2(1,:),Abot(1,:))
numel=10;
datamat2 = [[Xresc1(1,1:5:end)'; Xresc2(1,1:5:end)'], [tmpA(1,1:5:end)'; flip(Abot(1,1:5:end))']];

for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),tmpA(i,:))
    hold on
    plot(Xresc2(i,:),Abot(i,:))
    datamat2 = [datamat2, [[Xresc1(i,1:5:end)'; Xresc2(i,1:5:end)'], [tmpA(i,1:5:end)'; flip(Abot(i,1:5:end))']]];
end
figure; 
plot(Xresc1(1,:),[A0top(1,1); A0top(1,:)],'--')
hold on 
plot(Xresc2(1,:),Abot(1,:))

csvwrite('ADiscreteTimestepsgamma20.csv',datamat2); 

% PLOTTING U
figure; 
plot(Xresc1(1,:),ufull(1,:))
hold on 
plot(Xresc2(1,:),Abot(1,:))
numel=10;
datamat3 = [[Xresc1(1,1:10:end)'; Xresc2(1,1:10:end)'], [ufull(1,1:10:end)'; flip(Abot(1,1:10:end))']];

for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),ufull(i,:))
    hold on
    plot(Xresc2(i,:),Abot(i,:))
    datamat3 = [datamat3, [[Xresc1(i,1:10:end)'; Xresc2(i,1:10:end)'], [ufull(i,1:10:end)'; flip(Abot(i,1:10:end))']]];
end

csvwrite('uDiscreteTimestepsgamma20.csv',datamat3); 

return

% Round 2 
A0  = A(end,:);
%u0 = u(end,:);
th0 = th(end,:);
th0bot = flip( phi(end,:) ); 

y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = lam(end,:);  

y0(3*K+1:4*K) = th0bot;
% 

tic
[t,y] = ode15s(@coupledPde,tout,y0); 
toc
A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

lam   = y(:,2*K+1:3*K);

phi = y(:,3*K+1:4*K); 

u = zeros(size(A));
u0 = usolutionNoFB(A0',th0');  % Do I even need this guy? 
u(1,:) = u0; 
for i=2:N
    u(i,:) = usolution(A(i,:)',th(i,:)',lam(i,:)', L);   
    
end


% Abot  = y(:,4*K+1:5*K); % This should just be zeros


% We try again but with the achieved steady state from the previous method 

%% FOR NOW, IGNORE HOW WE OBTAIN THE SOLUTION TO U 
%% Solve for u at next time step (laplacian)
%     temp = 3*mu([0;0;th(:,i)]);     % add ghost node to th
%     tiph = ( temp(1:end-1) + temp(2:end) ) / 2;

%thfullT = [ zeros(N,1), ...
    %   thtop          ];
   
thtop = [ zeros(N,1), ...
       th         ];
  
Atop  = [ D*ones(N,1), ...
       A           ];
   
tmpA =  (Atop + [Atop(:,2:end), ones(N,1)])/2; % extract A at the edges 
Abot = ones(N,K); 
Afull =    [tmpA, Abot];  
ufull  = [ u , ...
       uf.*ones(N,1) ];
   
dx = L/K;

x = (0:dx:L)';
[X, T1] = meshgrid(x,t);
Xresc1 = [lam(:,1), lam].*X; 
[X, T2] = meshgrid((dx:dx:L)',t);
Xresc2 = -X.*(L-lam)+L;
set(gca,'TickLabelInterpreter','latex','fontsize',13)

% figure;
% surf(T1,Xresc1,thtop,'LineStyle','none')
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)


% figure; 
% surf(T2,Xresc2,phi,'LineStyle','none')
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% 
% figure; 
% surf(T1,Xresc1,thtop,'LineStyle','none')
% hold on 
% surf(T2,Xresc2,phi,'LineStyle','none')
% xlabel('$t$','Interpreter','latex')
% ylabel('$x$', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)

numel=20;
datamat = [[Xresc1(1,:)'; Xresc2(1,:)'], [thtop(1,:)'; phi(1,:)']];

figure; 
for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),thtop(i,:))
    hold on
    plot(Xresc2(i,:),phi(i,:)) 
    pause
    datamat = [datamat, [[Xresc1(i,:)';flip(Xresc2(i,:))'], [thtop(i,:)'; flip(phi(i,:))']]];
end
%csvwrite('ThetaDiscreteTimesteps.csv',datamat); 



numel=10;
datamat = [[Xresc1(1,:)'; flip(Xresc2(1,:))'], [tmpA(1,:)'; flip(Abot(1,:))']];
figure; 
for i = N/numel:(N/numel):N
    plot(Xresc1(i,:),tmpA(i,:))
    hold on
    plot(Xresc2(i,:),Abot(i,:))
    datamat = [datamat, [[Xresc1(i,:)';flip(Xresc2(i,:))'], [tmpA(i,:)'; flip(Abot(i,:))']]];
end
% 
% %csvwrite('ADiscreteTimesteps.csv',datamat); 
% 
% numel=10;
% datamat = [[Xresc1(1,:)'; flip(Xresc2(1,:))'], [ufull(1,:)'; flip(Abot(1,:))']];
% figure; 
% for i = N/numel:(N/numel):N
%     plot(Xresc1(i,:),ufull(i,:))
%     hold on
%     plot(Xresc2(i,:),Abot(i,:))
%     datamat = [datamat, [[Xresc1(i,:)';flip(Xresc2(i,:))'], [ufull(i,:)'; flip(Abot(i,:))']]];
% end
% csvwrite('uDiscreteTimesteps.csv',datamat);
% 
% this is a file that has two columns per timestep, where the first column
% is the x values and the second one the theta values. In the x and theta
% values, we have appended both sides onto one vector, so Xcol = Xresc1,
% Xresc2, and similarly for theta 



return 
figure;
% surf(t,x.*[lam(1,1); lam(:,end)],Atop','LineStyle','none')
% xlabel('$t$')
% ylabel('$X$')
% title('$A$ with MOL','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% figure;
% surf(t,x,thfull','LineStyle','none')
% xlabel('$t$')
% ylabel('$x$')
% title('$\theta$ with MOL','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)

figure;
surf(t,x,ufull','LineStyle','none')
xlabel('$t$')
ylabel('$X$')
title('$u$ with MOL','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

plot(t,lam(:,end))
xlabel('$t$','Interpreter','latex')
ylabel('$\lambda(t)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

datamat = [[t(1,end); t(10:10:end)],[lam(1,end); lam(10:10:end,end)]];

csvwrite('lambdavst.csv',datamat);

figure;
%surf(T1,Xresc1,Atop,'LineStyle','none')
surf(T3,Xresc3,Afull,'LineStyle','none')

% hold on
% Abot = ones(N,K); 
%surf(T2,Xresc2,Abot,'LineStyle','none')
xlabel('$t$','Interpreter','latex')
ylabel('$x$', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)



