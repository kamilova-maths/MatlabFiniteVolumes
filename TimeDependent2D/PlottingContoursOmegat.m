
xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

% PLOTTING THETA


% % Save data to file
figure;
title('Temperature')
tvector = omega*t*ones(1,2*K); 
% Temperature doesn't look quite there? A tiny cheat: (Shift the whole
% thing by 2*pi)
tvector = tvector - 2*pi*ones(size(tvector));
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,temp,20,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
xlim([0 2*pi])
if sav==1
    axis off
    print(gcf, '-dpng', '-r600', '-painters', 'Temperature.png')
end
% return 

% PLOTTING A
figure; 
tvector = omega*t*ones(1,2*K); 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,Acel, 20,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
xlim([0 2*pi])
if sav==1
    axis off
    print(gcf, '-dpng', '-r600', '-painters', 'A.png')
end
% PLOTTING U
figure; 
tvectorint = omega*t*ones(1,2*K+1);
xvectorint =  [lam*xint',lam + (L-lam)*xint(2:end)'];
contourf(tvectorint, xvectorint,uint,20,'LineColor', 'none')
if max(max(uint))>5
    maxu = 5
else
    maxu = max(max(uint))
end
caxis([min(min(uint)) maxu])
ax = gca;
ax.YDir = 'reverse';
xlim([0 2*pi])
if sav==1
    axis off
    print(gcf, '-dpng', '-r600', '-painters', 'Velocity.png')
    csvwrite('lam.csv',[omega*t, lam]);
   % csvwrite('P.csv',[t,P0t(t)]);
end
% PLOTTING lambda
figure; 
plot(omega*t, lam);
xlim([0 2*pi])
hold on
if P0tval==0
    plot(omega*t,P);  
    xlim([0 2*pi])
    if sav==1
        csvwrite('P.csv',[omega*t, P])
    end
elseif P0tval==1
    plot(omega*t,P0t(t));
    xlim([0 2*pi])
    if sav==1
        csvwrite('P0t.csv',[omega*t, P0t(t)])
    end
end
title('lambda and P')
xlabel('t')


if uftval ==1
    figure;
    plot(t,uft(t));
    xlim([0 2*pi])
    title('uf')
xlabel('t')
    if sav==1
        csvwrite('uft.csv',[t,uft(t)]);
    end
end




%code to save data 
