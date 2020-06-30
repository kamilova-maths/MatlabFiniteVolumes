

if st == 1 
    % Change filename to match what we want to import 
    %data = csvread('SSG23P0uin1oD.csv');
    x0cel = data(1:2*K);
    
    Asteady = data(2*K+1:4*K); 

    tempsteady=  data(4*K+1:6*K);
    
    x0int = data(6*K+1:8*K+1); 
    usteady = data(8*K+2:10*K+2); 
    lam0steady = data(10*K+3); 

end


% PLOTTING THETA

figure;
numel = 10; 
Kindices = [1, 5:5:2*K]; 
xvector1 = [xcel*lam(1);lam(1) + xcel*(L-lam(1))];
plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], temp(1,:)', '--');
thetadata = [xvector1(Kindices), temp(1,Kindices)'];
hold on

for i = floor(N/numel):floor((N/numel)):N
    xvector = [xcel*lam(i);lam(i) + xcel*(L-lam(i))];
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], temp(i,:)');
    thetadata = [thetadata, [xvector(Kindices), temp(i,Kindices)']];
end
%hold on 
if st == 1 
    plot(x0cel,tempsteady,'--')
    thetadata = [thetadata, [x0cel(Kindices), tempsteady(Kindices)]];
end

set(gca,'TickLabelInterpreter','latex','fontsize',13)
if dat == 1
csvwrite('ThetaDiscreteTimesteps.csv',thetadata); 
end
% 
% thsteadytop = temp(N,1:K);
% phisteady = temp(N,K+1:end);


title('Temperature')
% Save data to file

% return 
% IF you want to see the other solutions, remove the return. 
% PLOTTING A
figure; 
Adata = [xvector1(Kindices), Acel(1,Kindices)'];
plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], Acel(1,:)', '--');
	hold on
for i = floor(N/numel):floor((N/numel)):N
    xvector = [xcel*lam(i);lam(i) + xcel*(L-lam(i))];
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], Acel(i,:)');
    Adata = [Adata, [xvector(Kindices), Acel(i,Kindices)']];
end


if st == 1
    hold on
    plot(x0cel,Asteady,'--')
    Adata = [Adata, [x0cel(Kindices), Asteady(Kindices)]]; 
end

if dat == 1
csvwrite('ADiscreteTimesteps.csv',Adata); 
end
%Asteadytop = Asteady(1:300);
title('Area')

% Save data to file
%SS = [Asteadytop';thsteadytop';phisteady';lam(end,end)];
%csvwrite('SteadyStateK300.csv', SS); 

% PLOTTING U
figure; 
plot([xint*lam(1);lam(1) + xint(2:end)*(L-lam(1))], uint(1,:)', '--');
udata = [xvector1(Kindices), uint(1,Kindices)'];
hold on
for i = floor(N/numel):floor((N/numel)):N
    xvector = [xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))];
	plot([xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))], uint(i,:)');
    udata = [udata, [xvector(Kindices), uint(i,Kindices)']];
end


if st ==1
    hold on 
    plot(x0int,usteady,'--')
    udata = [udata, [x0int(Kindices), usteady(Kindices)]];
end

if dat == 1
    csvwrite('uDiscreteTimesteps.csv',udata); 
end

title('Velocity')
% Save data to file


% PLOTTING lambda
figure; 
plot(t, lam);
if dat == 1
csvwrite('lam.csv',[t, lam]);
%csvwrite('lamG23.csv',[[t(1), lam(1)]; [t(5:5:end), lam(5:5:end)]; [t(end), lam0steady]]);
end
title('lambda')
xlabel('t')

return

% PLOTTING Pt
figure; 
plot(t, P);
if dat == 1
csvwrite('P.csv',[t,P]);
%csvwrite('PG23.csv',[[t(1), P(1)]; [t(5:5:end), P(5:5:end)]; [t(end), P0]]);
end
title('P')
xlabel('t')

return

figure;
plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))],temp(1,:),'--')
hold on
plot([xcel*lam(N);lam(N) + xcel*(L-lam(N))],temp(N,:))

data = [[xcel*lam(1);lam(1) + xcel*(L-lam(1))], temp(1,:)', ...
    [xcel*lam(N);lam(N) + xcel*(L-lam(N))], temp(N,:)']; 

csvwrite('ThetaSteadyComparisonP0tG30.csv',data);

figure;
plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))],Acel(1,:),'--')
hold on
plot([xcel*lam(N);lam(N) + xcel*(L-lam(N))],Acel(N,:))

data = [[xcel*lam(1);lam(1) + xcel*(L-lam(1))], Acel(1,:)', ...
    [xcel*lam(N);lam(N) + xcel*(L-lam(N))], Acel(N,:)']; 

csvwrite('ASteadyComparisonP0tG30.csv',data);

figure;
plot([xint*lam(1);lam(1) + xint(2:end)*(L-lam(1))],uint(1,:),'--')

hold on
plot([xint*lam(N);lam(N) + xint(2:end)*(L-lam(N))],uint(N,:))

data = [[xint*lam(1);lam(1) + xint(2:end)*(L-lam(1))], uint(1,:)', ...
   [xint*lam(N);lam(N) + xint(2:end)*(L-lam(N))], uint(N,:)']; 

csvwrite('uSteadyComparisonP0tG30.csv',data);

%% NOTE: I don't understand the video version ... but maybe it would be nice for a presentation? 
% video version
figure('units','normalized','outerposition',[0 0 0.25 1])
for i = N/numel:(N/numel):N
	subplot(3,1,1)
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], temp(i,:)'), axis([0 L 0 (max(max(temp))+0.1)]);
	hold on
	plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], temp(1,:)', '--'), axis([0 L 0 (max(max(temp))+0.1)]);
	plot([lam(i),lam(i)], [0,max(max(temp))+0.1]);
	hold off
	title(strcat('T=',num2str(i*T/N)))
	ylabel('Temperature')
	subplot(3,1,2)
	plot([xcel*lam(i);lam(i) + xcel*(L-lam(i))], Acel(i,:)'), axis([0 L 0 (max(max(Acel))+0.1)]);
	hold on
	plot([xcel*lam(1);lam(1) + xcel*(L-lam(1))], Acel(1,:)', '--'), axis([0 L 0 (max(max(Acel))+0.1)]);
	hold off
	ylabel('Area')
  subplot(3,1,3)
	plot([xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))], uint(i,:)'), axis([0 L 0 (max(max(uint))+0.1)]);
	hold on
	plot([xint*lam(i);lam(i) + xint(2:end)*(L-lam(i))], uint(i,:)', '--'), axis([0 L 0 (max(max(uint))+0.1)]);
	hold off
	ylabel('Velocity')
	pause(T/N)
    %pause(1)
end


