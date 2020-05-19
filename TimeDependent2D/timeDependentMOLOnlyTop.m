% Clear previous files
clear all 
close all 
clc

% Parameters shared with other routines (I'm not very sure about this,but
% let's give it a try)
global Pe Bi tha N K gamma P0 St T L D uf 

% Define parameters
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
%St=(rho*g*L^2)/(uc*mu0);
St = 21; 
P0 = (10000*L)/((R1^2)*uc*mu0);

DeltaT = (Qc*L)/(rho*c*uc);
Bi= ((L^2)*h)/(k*R1); 
tha = 0.005; 
D = (R0^2)/(R1^2); 

gamma = 0; 

%This is the area of the clamps, taken from Temperature profiles ... 
x1 = 5/7;
x2 = 6.5/7;
Q = 1;
eps = 1e-2;

% Calculating the initial conditions as a solution of the steady state
% problem 
N=1000; K=2000;
% end of the domain
T = 1; L= 1 ;
dx = L/K;

% We add the heaviside with H=1, and we remove it with H=0. 
H=0;

% Plots for steady state - 1 , no plots for steady state - 0
plt = 0;

% We try with tha=0
%[P, A0, J0, th0, ~] = InitialConditionsSteady(K,gamma,Q,x1,x2,eps,St,tha,Bi,Pe,P0,R0,L,H,plt);

% Analytic solution for A and P (and u?) 0 for the steady state, with
% constant mu (no temperature) 
lamex = @(St) atan(sqrt(P0^2 + 6*St - 6*D*St)/sqrt(6*St*D-P0^2)).*6./(sqrt(6*St*D-P0^2)) - ...
    6*atan(P0./sqrt(6*St*D-P0^2))./(sqrt(6*St*D-P0^2)); 

newSt = fsolve(@(St) lamex(St)-1, 21);
St = newSt;

x = linspace(0,L,K);
fac = 6*St*D-P0^2;
P0bar =  6*atan(P0./sqrt(fac))./(sqrt(fac)); 
Aex = (fac/(6*St)).*tan(sqrt(fac).*(x+P0bar)./6).^2 - P0.^2./(6*St)+D ; 


lam0 = lamex(St); 
% Note that Aex is A evaluated at the nodes. I want to evaluate it at the
% cells, so I have to perform an averaging, of the form 
Aex = [D; Aex'];
A0 = (Aex(1:end-1) + Aex(2:end))./2; 

% Find lam0, then resize both sides with an interpolation (only necessary
% to do this for the steady state, the conditions on the rest are much
% nicer because of our rescalings 
%Initial conditions given by A0, th0
% A0=A0';
% th0=th0'; 
% lam0 = 1; 

th0 = zeros(size(A0)); 
uf = 1/A0(end); 
u0 = usolution(A0,th0,lam0); 


y0(1:K) = A0;
y0(1+K:2*K) = lam0.*ones(size(A0));

% Independent variable for ODE integration 
tout = linspace(0,T,N);

% ODE integration 
reltol = 1.0e-04; abstol = 1.0e-04;
options = odeset('RelTol',reltol,'AbsTol',abstol);
[t,y] = ode15s(@coupledPdeNoTemp,tout,y0);
%[t,y] = ode15s(@coupledPde,tout,y0); 

A  = y(:,1:K);
lam   = y(:,1+K:2*K);

%% Solve for u at next time step (laplacian)
% We construct u from A 

umat = zeros(size(A));
umat(1,:) = u0; 
for i=2:N
    umat(i,:) = usolution(A(i,:)',zeros(size(A(i,:)))',lam(i,:)'); 
end

Afull  = [ D*ones(N,1), ...
       A          ];   
   
% We rescale it as 

   
ufull  = [ umat , ...
       uf.*ones(N,1) ];
   
dx = 1/K;

x = (0:dx:L)';
A0 = [D; A0]; 
figure; 
plot(x,A0,'--')
hold on
plot(x,Afull(100:100:end,:))
set(gca,'TickLabelInterpreter','latex','fontsize',13)
xlabel('$x$','Interpreter','latex')
ylabel('$A$','Interpreter','latex')
% A0 = [D; A0]; 
figure; 
u0 = [u0; uf];
plot(x,u0,'--')
hold on
plot(x,ufull(100:100:end,:))
set(gca,'TickLabelInterpreter','latex','fontsize',13)
xlabel('$x$','Interpreter','latex')
ylabel('$u$','Interpreter','latex')
% figure;
% surf(t,x,Afull','LineStyle','none')
% xlabel('$t$')
% ylabel('$X$')
% title('$A$ with MOL')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% figure;
% contourf(t,x,Afull','LineStyle','none')
% colorbar 
% xlabel('$t$')
% ylabel('$x$')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 

figure;
plot(t,lam(:,end))
%title('True $\lambda = $')
xlabel('$t$','Interpreter','latex')
ylim([0.5 1])
ylabel('$\lambda(t)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

return
valuesmatrix=[x(1:10:end),A0(1:10:end),Afull(100:100:end,1:10:end)'];
csvwrite('Awithlambda1constantmu.csv',valuesmatrix)

valuesmatrix=[t(1:10:end),lam(1:10:end,end)] ; 
csvwrite('lambdaconstantmu.csv', valuesmatrix)

valuesmatrix=[x(1:10:end),u0(1:10:end),ufull(100:100:end,1:10:end)'];
csvwrite('uwithlambda1constantmu.csv',valuesmatrix)

% figure;
% surf(t,x,thfull')
% xlabel('$t$')
% ylabel('$x$')
% title('$\theta$ with MOL')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% 
% figure;
% surf(t,x,ufull')
% xlabel('$t$')
% ylabel('$x$')
% title('$u$ with MOL')
% set(gca,'TickLabelInterpreter','latex','fontsize',13)

