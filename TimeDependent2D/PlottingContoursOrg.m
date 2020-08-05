
xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
%indx = N;
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
%tvector = t*ones(1,2*K); 
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
if max(max(uint))>5
    maxu = 5
else
    maxu = max(max(uint))
end
caxis([min(min(uint)) maxu])
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
    plot(t, P);
    if sav==1
        csvwrite('P.csv',[t, P])
    end
elseif P0tval==1
    plot(t,P0t(t));
    if sav==1
        csvwrite('P0t.csv',[t, P0t(t)])
    end
end
title('lambda and P')
xlabel('t')

if uftval ==1
    figure;
    plot(t,uft(t));
    title('uf')
    xlabel('t')
    if sav==1
        csvwrite('uft.csv',[omega*t,uft(t)]);
    end
end


return

%code to calculate average P0


indx1 = find(t==first(step),1);



space = length(first(step:end));
space = space - 1; 
indx2 = find(t==first(step+space),1);
tbetween = t(indx1:indx2);

[tbetween,index] = unique(tbetween);
Pinterp = interp1(tbetween,P(index),linspace(t(indx1),t(indx2),length(tbetween)));
mean(Pinterp)
Paverage = sum(P(indx1:indx2))/(indx2-indx1)

