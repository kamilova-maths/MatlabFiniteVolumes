%% START FROM HERE WHEN YOU HAVE ALREADY CALCULATED STEADY STATE 
%close all
y0(1:K) = A0;
y0(1+K:2*K) = A0.*th0;
y0(2*K+1:3*K) = th0bot;  

y0(3*K+1) = lam0;

% Independent variable for ODE integration 
tout = linspace(0,T,N);

% ODE integration 
reltol = 1.0e-04; abstol = 1.0e-04;
options = odeset('RelTol',reltol,'AbsTol',abstol);
%[t,y] = ode15s(@coupledPdeNoTemp,tout,y0);
tic
[t,y] = ode15s(@coupledPde,tout,y0); 
toc
%pause
A  = y(:,1:K); % This is A from X=0 to X=1 (this is, 0<x<lambda)

th = y(:,K+1:2*K)./A;

phi   = y(:,2*K+1:3*K);

lam = y(:,3*K+1); 

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
top = max(max(abs(th(1,:)-th(end,:))))
bottom = max(max(abs(phi(1,:)-phi(end,:))))