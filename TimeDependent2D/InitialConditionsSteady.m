function [P, A, J, Th, x,lam] = InitialConditionsSteady(eps,H,plt)
% DATE:     2020 
% DESCR:    [P, A, J, Th, x,lam] = InitialConditionsSteady(eps,H,plt)
%
%           Simple code that uses bvp4c to calculate the steady state 
%           solution to the coupled flow and temperature problem. It is
%           necessary to run ParametersDefinition before running this code.
% INPUT: 
%           eps: (scalar) small value that determines the transition for the
%           Heaviside function used for A solution. An acceptable value is
%           O(10^-4) or smaller, normally. 
%           H:  (scalar) indicator for use of Heaviside. If H == 1, we use the
%           version of the problem with a Heaviside function. Otherwise,
%           A is allowed to increase past 1.
%           plt: (scalar) indicator for plotting. plt==1 enables plotting.
%           Plotting is disabled otherwise. 
%          
% OUTPUT:   P:  (row vector) Pressure, corresponding to -3 mu(theta) du/dx. 
%           A:  (row vector) Area, changes from D to 1, and then remains
%           constant (when H ==1 ).
%           J:  (row vector) Heat flux, comes out from rewriting the
%           equation for theta to make it into two first order ODEs. We
%           only have a condition at the end of the domain for J.
%           Th: (row vector) Temperature. 
%           All of P, A, J, and Th are of the same size. The length of
%           these row vectors is automatically determined by bvp4c, as it 
%           uses adapts the mesh according to the computed gradients. 
%           x : (row vector) Contains the x values for the computed
%           solutions of the ODEs. 
%           lam: (scalar) Position where A reaches 1 for the first time. 
% ADDITIONAL COMMENTS: 
%           There is a twin version of this code in the steady
%           state folder. Leave both there. I think it's worth having a 
%           functional example in the Time Dependent part as well. 
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           bvpinit, bvp4c: to solve the ODEs
%           odefun: Contains the equations to solve, for H == 1
%           odefunNH: Contains the equations to solve, for H \neq 1
%           bcfun: Contains the boundary conditions for the ODEs. 


% Invoking global variables. 
global L

% Initial guess for the solution. These are just constants, more or less
% corresponding to the average values. These are not necessarily the
% boundary conditions. 
Tin = 0;
Pin = 1;
Ain = 0.25;

% Initial solution using the initial guess. 

solinit = bvpinit(linspace(0,L,100),@(x) [Pin; Ain; 0; Tin]);

% We set the tolerances for the solution. 
opts = bvpset('RelTol',1e-3,'AbsTol',1e-6);

% We choose to either use or not use the Heaviside, depending upon input
% parameter H . This is where the solution is computed.
if H==1
    sol = bvp4c(@(x,y) odefun(x,y,eps),@(ya,yb) bcfun(ya,yb,Tin),solinit,opts);
else
	sol = bvp5c(@(x,y) odefunNH(x,y),@(ya,yb) bcfun(ya,yb,Tin),solinit,opts);
end

% We extract the solutions. Each one of these is a row vector of length(x).
% 
P  = sol.y(1,:);
A  = sol.y(2,:);
J  = sol.y(3,:); 
Th = sol.y(4,:); 
x = sol.x; 


% If plotting is selected, we proceed to this if statement. 
if plt==1
% We take this long and convoluted way to plot the green shaded region where 
% Q is activated.     
xplt = linspace(-0.5,1.1,1000); 
X = [xplt,fliplr(xplt)]; 
y1 = x1*ones(size(xplt)); 
y2 = x2*ones(size(xplt));
Y = [y1, fliplr(y2)]; 

fill(X,Y,[204/255, 255/255, 204/255], 'LineStyle','none','HandleVisibility','off');
hold on 
figure(1)
plot(sol.y(2:end,:),sol.x); drawnow; shg;
set(gca,'YDir','reverse')
ylabel('$x$','Interpreter','latex', 'Fontsize', 15)
set(gca,'YDir','reverse','TickLabelInterpreter','latex','fontsize',15)
hold on 
% We take another convoluted way to plot a dashed line where lambda is. 
i = find(sol.y(2,:)>0.9999,1);
lam=sol.x(i);
plot(x,sol.x(i)*ones(size(x)),'--k'); drawnow; shg;
leg1 = legend('$A(x)$', '$J(x)$','$\theta(x)$','$x=\lambda$');
set(leg1,'Interpreter','latex','fontsize',14)
xlim([0, 1.1])
hold off
% 

end

end