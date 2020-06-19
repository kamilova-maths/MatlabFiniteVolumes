
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
csvwrite('lam.csv',[[t(1), lam(1)]; [t(5:5:end), lam(5:5:end)]]);
hold on
plot(t,P0t(t)); 
csvwrite('P0t.csv',[[t(1), P0t(t(1))]; [t(5:5:end), P0t(t(5:5:end))]]); 
%plot([t(1),t(end)], [lam(1),lam(1)]);
title('lambda')
xlabel('t')

%code to save data 
