function SolvingSteadyState()
% 3 October 2018

% Modified version of Ian Hewitt's original code. Matches transfer of
% status work. 

% Realistic values were included for the data that is available from
% papers. 

close all;

%Data 

rho= 1.8*10^3; %Bergstrom ; 
g = 10; 
c= 900; % Fitt and Howell
L=7; % Temperature Profiles in Soderberg Electrodes
u = 10^-5; %Bergstrom approximation
R0=0.5;
R1=1; 
k = 3; 
Qc=15000; %Taken very vaguely from Temperature profiles in Soderberg electrodes. 
mu0=10^10; % given by Bjornar at a reference temperature
T_a=343; 
h = 7; 

%Defining non-dimensional parameters
% Peclet number
p.Pe = (rho*c*u*L)/(k);
p.epsilon=R1/L;
p.St=(rho*g*L^2)/(u*mu0);
p.P0 = (10000*L)/((R1^2)*u*mu0);

p.DeltaT = (Qc*L)/(rho*c*u);

p.Bi= ((L^2)*h)/(k*R1); 


p.alpha = 50; 
p.Ta = 0;
p.Tin = 0;
p.Pin = 0;

%p.Pin = 0; %produces nicer plots
p.Ain = 0.5;

%This is the area of the clamps, taken from Temperature profiles ... 
p.x1 = 5/7;
p.x2 = 6.5/7;
p.Q = 1;
p.eps = 1e-2;
p1 = p;

% % initial solve using easy parameters

solinit = bvpinit(linspace(0,1,100),@(x) [p.Pin; p.Ain; 0; p.Tin]);
opts = bvpset('RelTol',1e-4,'AbsTol',1e-4);
sol = bvp4c(@odefun,@bcfun,solinit,opts);


%plot(xlim, [1 1]*p.x1, '--')
%plot(xlim, [1 1]*p.x2, '--')
figure(1)
x = linspace(-0.5,1.1,1000); 
X = [x,fliplr(x)]; 
y1 = p.x1*ones(size(x)); 
y2 = p.x2*ones(size(x));
Y = [y1, fliplr(y2)]; 

fill(X,Y,[204/255, 255/255, 204/255], 'LineStyle','none','HandleVisibility','off');
hold on 

plot(sol.y(2:end,:),sol.x); drawnow; shg;
set(gca,'YDir','reverse')
ylabel('$x$','Interpreter','latex', 'Fontsize', 15)
set(gca,'YDir','reverse','TickLabelInterpreter','latex','fontsize',15)
hold on 
i = find(sol.y(2,:)>0.9999,1);
sol.x(i)
plot(x,sol.x(i)*ones(size(x)),'--k'); drawnow; shg;
leg1 = legend('$A(x)$', '$J(x)$','$\theta(x)$','$x=\lambda$');
set(leg1,'Interpreter','latex','fontsize',14)
xlim([-0.5, 1.1])
hold off

figure(2)
x = linspace(-2,6,1000); 
X = [x,fliplr(x)]; 
y1 = p.x1*ones(size(x)); 
y2 = p.x2*ones(size(x));
Y = [y1, fliplr(y2)]; 

% laplacian for theta
Dx2 = spdiags( ones(N,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], N, N );
Dx2(N,N-1) = 2/dx^2;


% solve for u0
fill(X,Y,[204/255, 255/255, 204/255], 'LineStyle','none','HandleVisibility','off');
hold on 
plot(sol.y(1,:),sol.x); drawnow; shg;
set(gca,'YDir','reverse')
ylabel('$x$','Interpreter','latex', 'Fontsize', 15)
set(gca,'YDir','reverse','TickLabelInterpreter','latex','fontsize',15)
hold on 
i = find(sol.y(2,:)>0.9999,1);
sol.x(i)
plot(x,sol.x(i)*ones(size(x)),'--k'); drawnow; shg;
leg1 = legend('$P(x)$','$x=\lambda$');
set(leg1,'Interpreter','latex','fontsize',14)
hold off
    
function dydx = odefun(x,y)
P = y(1,:);
A = min(y(2,:),1);
T = y(4,:);
J = y(3,:);
mu = exp(-p.alpha*T);  
% viscosity
Q = p.Q*(x>p.x1).*(x<p.x2);    % heat source
Ta = p.Ta;                     % air temperature
%He = (A<1);                     % Heaviside function
% He = 1/2*(1+tanh((1-A)/p.eps)); % regularised Heaviside function
He = max(1-exp(-(1-A)/p.eps),0);% regularised Heaviside function
dydx = [ p. St*A;
        (A.*P./(3*mu)).*He;
        p.Pe*((J./A) - A.*Q) - 2*p.Bi*A.^(1/2).*(Ta-T);
        J./A
        ];
end

function res = bcfun(ya,yb)
res = [ ya(1,:)-p.P0;
        ya(2,:)-R0^2; 
        yb(3,:);
        ya(4,:)-p.Tin 
        ];
end
% function yinit = odefuninit(sol)
%        yinit = sol.y; 
% end 
end