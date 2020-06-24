
xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

% PLOTTING THETA


% % Save data to file
figure;
title('Temperature')
tvector = t*ones(1,2*K); 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,temp,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
axis off
print(gcf, '-dpng', '-r600', '-painters', 'Temperature.png')

% return 

% PLOTTING A
figure; 
tvector = t*ones(1,2*K); 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,Acel,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
axis off
print(gcf, '-dpng', '-r600', '-painters', 'A.png')

% PLOTTING U
figure; 


tvectorint = t*ones(1,2*K+1);
xvectorint =  [lam*xint',lam + (L-lam)*xint(2:end)'];
contourf(tvectorint, xvectorint,uint,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
axis off
print(gcf, '-dpng', '-r600', '-painters', 'Velocity.png')
% PLOTTING lambda
figure; 
plot(t, lam);
csvwrite('lamG23.csv',[[t(1), lam(1)]; [t(5:5:end), lam(5:5:end)]; [t(end), lam0]]);
hold on
plot(t,P); 
csvwrite('PG23.csv',[[t(1), P(1)]; [t(5:5:end), P(5:5:end)]; [P(end), P0]]);
%plot([t(1),t(end)], [lam(1),lam(1)]);
title('lambda and P')
xlabel('t')

%code to save data 
