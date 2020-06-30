function [P, A, J, Th, x] = InitialConditionsSteady(eps,H,plt)
% 3 October 2018

% Modified version of Ian Hewitt's original code. Matches transfer of
% status work. 

% Realistic values were included for the data that is available from
% papers. 
% We increase the size of the discretisation to avoid extrapolation (which
% I think might cause me issues in future)
global x1 x2 L
close all;

%Data
% Indicator for plotting. 0- no plots, 1- plots

Tin = 0;
Pin = 1;

Ain = 0.25;

% % initial solve using easy parameters

solinit = bvpinit(linspace(0,L,100),@(x) [Pin; Ain; 0; Tin]);
opts = bvpset('RelTol',1e-3,'AbsTol',1e-6);

if H==1
    sol = bvp4c(@(x,y) odefun(x,y,eps),@(ya,yb) bcfun(ya,yb,Tin),solinit,opts);
else
	sol = bvp5c(@(x,y) odefunNH(x,y,eps),@(ya,yb) bcfun(ya,yb,Tin),solinit,opts);
end

P  = sol.y(1,:);
A  = sol.y(2,:);
J  = sol.y(3,:); 
Th = sol.y(4,:); 
x = sol.x; 

%Ainterp=interp1(sol.x,sol.y(2,:),linspace(0,L,n),'spline');

% option #1 - taking the values evaluated at the point


% option #2 - taking the values evaluated at the cell
%uf=1


if plt==1
x = linspace(-0.5,1.1,1000); 
X = [x,fliplr(x)]; 
y1 = x1*ones(size(x)); 
y2 = x2*ones(size(x));
Y = [y1, fliplr(y2)]; 

fill(X,Y,[204/255, 255/255, 204/255], 'LineStyle','none','HandleVisibility','off');
hold on 
figure(1)
plot(sol.y(2:end,:),sol.x); drawnow; shg;
set(gca,'YDir','reverse')
ylabel('$x$','Interpreter','latex', 'Fontsize', 15)
set(gca,'YDir','reverse','TickLabelInterpreter','latex','fontsize',15)
hold on 
i = find(sol.y(2,:)>0.9999,1);
sol.x(i);
plot(x,sol.x(i)*ones(size(x)),'--k'); drawnow; shg;
leg1 = legend('$A(x)$', '$J(x)$','$\theta(x)$','$x=\lambda$');
set(leg1,'Interpreter','latex','fontsize',14)
xlim([-0.5, 1.1])
hold off
% 
% figure(2)
% x = linspace(-2,6,1000); 
% X = [x,fliplr(x)]; 
% y1 = p.x1*ones(size(x)); 
% y2 = p.x2*ones(size(x));
% Y = [y1, fliplr(y2)]; 
% fill(X,Y,[204/255, 255/255, 204/255], 'LineStyle','none','HandleVisibility','off');
% hold on 
% plot(sol.y(1,:),sol.x); drawnow; shg;
% set(gca,'YDir','reverse')
% ylabel('$x$','Interpreter','latex', 'Fontsize', 15)
% set(gca,'YDir','reverse','TickLabelInterpreter','latex','fontsize',15)
% hold on 
% i = find(sol.y(2,:)>0.9999,1);
% sol.x(i)
% plot(x,sol.x(i)*ones(size(x)),'--k'); drawnow; shg;
% leg1 = legend('$P(x)$','$x=\lambda$');
% set(leg1,'Interpreter','latex','fontsize',14)
% hold off
end




% function yinit = odefuninit(sol)
%        yinit = sol.y; 
% end 
end