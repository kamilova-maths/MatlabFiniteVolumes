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
N=800; K=300;
% end of the domain
T = 5; L= 1 ;
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
u0 = 1./Aex'; 
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
uf = 1;
%u0 = 1./A0; 
%u0 = usolution(A0,th0,lam0,L,P0); 


% We now check the return to steady state 
y0(1:K) = (1- 1e-2 -D).*linspace(0,1,K)'+D; 

y0(1+K) = 0.8;

% Independent variable for ODE integration 
tout = linspace(0,T,N);

% ODE integration 
reltol = 1.0e-04; abstol = 1.0e-06;
options = odeset('RelTol',reltol,'AbsTol',abstol);
[t,y] = ode15s(@coupledPdeNoTemp,tout,y0);
%[t,y] = ode15s(@coupledPde,tout,y0); 

A  = y(:,1:K);
lam   = y(:,1+K);

%% Solve for u at next time step (laplacian)
% We construct u from A 

u = zeros(size(A));
%umat(1,:) = u0; 
for i=1:N
    u(i,:) = usolution(A(i,:)',zeros(size(A(i,:)))',lam(i),L,P0); 
end

Acel = [ A, ones(N,K)];		

%Afull  = [ D*ones(N,1), ...
%       A          ];   
   
% We rescale it as 

uint  = [ u , ...
		uf.*ones(N,K+1) ];   
   
dx = 1/K;

x = (0:dx:L)';
% A0 = [D; A0]; 
% figure; 
% plot(x,A0,'--')
% hold on
% plot(x,Afull(100:100:end,:))
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% xlabel('$x$','Interpreter','latex')
% ylabel('$A$','Interpreter','latex')

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
xvector1 = [xcel*lam(1);lam(1) + xcel*(L-lam(1))];
figure; 
numel = 10; 
Kindices = [1, 5:5:2*K]; 

Adata = [xvector1(Kindices), Acel(1,Kindices)'];
plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], Acel(1,:)', '--');
	hold on
for i = N/numel:(N/numel):N
    xvector = [xcel*lam(i);lam(i) + xcel*(L-lam(i))];
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], Acel(i,:)');
    Adata = [Adata, [xvector(Kindices), Acel(i,Kindices)']];
end
csvwrite('Awithlambda1constantmu.csv',Adata); 
hold on 
xsteadyint = linspace(0,lam0,K)';
xsteadycel = linspace(xsteadyint(2)/2,1-xsteadyint(2)/2,K)'; 
plot(xsteadycel,A0,'--')
csvwrite('ASteadywithlambda1constantmu.csv', [xsteadycel([1, 5:5:K]), A0([1, 5:5:K])])

% A0 = [D; A0]; 

% Plot it the right way, K
figure; 
plot([xint*lam(1);lam(1) + xint(2:end)*(L-lam(1))], uint(1,:)', '--');
udata = [xvector1(Kindices), uint(1,Kindices)'];
hold on
for i = N/numel:(N/numel):N
    xvector = [xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))];
	plot([xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))], uint(i,:)');
    udata = [udata, [xvector(Kindices), uint(i,Kindices)']];
end
plot(linspace(0,lam0,K),u0,'--')
csvwrite('uwithlambda1constantmu.csv',udata); 
csvwrite('uSteadywithlambda1constantmu.csv', [xsteadyint([1, 5:5:K]), u0([1, 5:5:K])])

% figure; 
% xu = linspace(0,L,K+1);
% %u0 = [u0; uf];
% plot(xu,ufull(1,:),'--')
% hold on
% plot(xu,ufull(100:100:end,:))
% set(gca,'TickLabelInterpreter','latex','fontsize',13)
% xlabel('$x$','Interpreter','latex')
% ylabel('$u$','Interpreter','latex')
% 
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
plot(t,lam)
%title('True $\lambda = $')
xlabel('$t$','Interpreter','latex')
ylim([0.5 1])
ylabel('$\lambda(t)$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',13)

csvwrite('lambdaconstantmu.csv', [t,lam])
return 
valuesmatrix=[x(1:10:end),Afull(1,1:10:end)',Afull(100:100:end,1:10:end)',A0(1:10:end)];
csvwrite('Awithlambda1constantmu.csv',valuesmatrix)



valuesmatrix=[xu(1:10:end)',ufull(1,1:10:end)',ufull(100:100:end,1:10:end)',[u0(1:10:end); uf]];
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

