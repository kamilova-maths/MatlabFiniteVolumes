
xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
%indx = N;
% PLOTTING THETA
patt = 7; 
rep  = 2; 
step = patt*rep*d;

list = patt*rep:floor(T/d);
dvalues = d*list; 
dvalues = dvalues-patt*rep*d; 
cylvec = [dvalues(1:end-1)', cyladd(patt*rep:end)] ; 
csvwrite('DStCyl.csv', cylvec); 
njvalues = [dvalues(1:end-1)',cylvec(:,2)./onecyl];
csvwrite('njvalues.csv',njvalues);
csvwrite('height.csv',[njvalues(:,1), c1dim*njvalues(:,2)])

% Save data to file
figure;
title('Temperature')
%tvector = t*ones(1,2*K); 
tvector = t*ones(1,2*K)-step; 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,temp,20,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
xlim([0 tvec(end)-step])
if sav==1
    axis off
    print(gcf, '-dpng', '-r600', '-painters', 'Temperature.png')
end


%PLOTTING A
figure; 
%tvector = t*ones(1,2*K)-step; 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
contourf(tvector, xvector,Acel, 20,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
xlim([0 tvec(end)-step])
if sav==1
    axis off
    print(gcf, '-dpng', '-r600', '-painters', 'A.png')
end
%PLOTTING U
figure; 
tvectorint = t*ones(1,2*K+1)-step;
xvectorint =  [lam*xint',lam + (L-lam)*xint(2:end)'];
contourf(tvectorint, xvectorint,uint,80,'LineColor', 'none')
xlim([0 tvec(end)-step])
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
    csvwrite('lam.csv',[t-step, lam]);
end
% PLOTTING lambda
figure; 
plot(t-step, lam);
xlim([0 t(end)-step])

hold on
if P0tval==0
    plot(t-step,P);
     xlim([0 t(end)-step])
    if sav==1
        csvwrite('P.csv',[t-step, P])
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
        csvwrite('uft.csv',[t,uft(t)]);
    end
end


return

%code to calculate average P0


%indx1 = find(t==first(step),1);

indx1 = find(t==step,1);
indx2 = find(tvec == patt*(rep+1)*d,1);

[tunique, index] = unique(tvec);

Pinterp = interp1(tunique,P(index),linspace(t(indx1),t(indx2),100));
mean(Pinterp)
Paverage = sum(P(indx1:indx2))/(indx2-indx1)

