
xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

% PLOTTING THETA


% % Save data to file
figure;
title('Temperature')
tvector = t*ones(1,2*K); 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,temp,20,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
if sav==1
    axis off
    print(gcf, '-dpng', '-r600', '-painters', 'Temperature.png')
end
% return 

% PLOTTING A
figure; 
tvector = t*ones(1,2*K); 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,Acel, 20,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
if sav==1
    axis off
    print(gcf, '-dpng', '-r600', '-painters', 'A.png')
end
% PLOTTING U
figure; 
tvectorint = t*ones(1,2*K+1);
xvectorint =  [lam*xint',lam + (L-lam)*xint(2:end)'];
contourf(tvectorint, xvectorint,uint,20,'LineColor', 'none')
%caxis([1 5])
ax = gca;
ax.YDir = 'reverse';
if sav==1
    axis off
    print(gcf, '-dpng', '-r600', '-painters', 'Velocity.png')
    csvwrite('lam.csv',[t, lam]);
   % csvwrite('P.csv',[t,P0t(t)]);
end
% PLOTTING lambda
figure; 
plot(t, lam);

hold on
if P0tval==0
    plot(t,P);  
    csvwrite('P.csv',[t, P])
else
    plot(t,P0t(t));
    csvwrite('P0t.csv',[t, P0t(t)])
end
%plot(t,P0t(t));
%plot([t(1),t(end)], [lam(1),lam(1)]);
title('lambda and P')
xlabel('t')

%code to save data 
