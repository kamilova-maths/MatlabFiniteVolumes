
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
set(gca,'TickLabelInterpreter','latex','fontsize',13)
hold on 
plot(linspace(0,lam0steady,K),th0steady,'--')
hold on 
plot(linspace(lam0steady,L,K),phi0steady,'--')


title('Temperature')
% % Save data to file
figure;
tvector = t*ones(1,2*K); 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,temp,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
axis off


% return 
% IF you want to see the other solutions, remove the return. 
% PLOTTING A
figure; 
tvector = t*ones(1,2*K); 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,Acel,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
axis off


% PLOTTING U
figure; 


tvectorint = t*ones(1,2*K+1);
xvectorint =  [lam*xint',lam + (L-lam)*xint(2:end)'];
contourf(tvectorint, xvectorint,uint,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
axis off

% Save data to file


% PLOTTING lambda
figure; 
plot(t, lam);
csvwrite('lam.csv',[[t(1), lam(1)]; [t(5:5:end), lam(5:5:end)]]);
hold on
plot(t,P0t(t)); 
csvwrite('P0t.csv',[[t(1), P0t(t(1))]; [t(5:5:end), P0t(t(5:5:end))]]); 
%plot([t(1),t(end)], [lam(1),lam(1)]);
title('lambda')
xlabel('t')

