function [Pinterp, Ainterp, Jinterp, Thinterp,uf] = InitialConditionsSteady(n,gamma,Q,x1,x2,eps,St,Ta,Bi,Pe,P0,R0,L,H,plt,lam)
% 3 October 2018

% Modified version of Ian Hewitt's original code. Matches transfer of
% status work. 

% Realistic values were included for the data that is available from
% papers. 
% We increase the size of the discretisation to avoid extrapolation (which
% I think might cause me issues in future)

close all;

%Data
% Indicator for plotting. 0- no plots, 1- plots

Tin = 0;
Pin = 0;

Ain = 0.2;

% % initial solve using easy parameters

solinit = bvpinit(linspace(0,L,100),@(x) [Pin; Ain; 0; Tin]);
opts = bvpset('RelTol',1e-6,'AbsTol',1e-7);

if H==1
    sol = bvp5c(@(x,y) odefun(x,y,gamma,Q,x1,x2,eps,St,Ta,Bi,Pe),@(ya,yb) bcfun(ya,yb,P0,R0,Tin),solinit,opts);
else
	sol = bvp5c(@(x,y) odefunNH(x,y,gamma,Q,x1,x2,eps,St,Ta,Bi,Pe),@(ya,yb) bcfun(ya,yb,P0,R0,Tin),solinit,opts);
end

%Pinterp  = sol.y(1,:);
%Ainterp  = sol.y(2,:);
%Thinterp = sol.y(4,:); 
Pinterp=interp1(sol.x,sol.y(1,:),linspace(0,L,n),'spline');
Ainterp=interp1(sol.x,sol.y(2,:),linspace(0,L,n+1),'pchip');
Thinterp = interp1(sol.x,sol.y(4,:),linspace(0,L,n+1),'pchip');
%Ainterp=interp1(sol.x,sol.y(2,:),linspace(0,L,n),'spline');

% option #1 - taking the values evaluated at the point
 uf=1/Ainterp(end);

Ainterp=(Ainterp(1:end-1)+Ainterp(2:end))/2;                % let's recover cell-averages for this
Thinterp=(Thinterp(1:end-1)+Thinterp(2:end))/2;

% option #2 - taking the values evaluated at the cell
%uf=1/Ainterp(end);

Jinterp=interp1(sol.x,sol.y(3,:),linspace(0,L,n));
%Thinterp=interp1(sol.x,sol.y(4,:),linspace(0,L,n),'spline');


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