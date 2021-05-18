% DATE:     2020 
% DESCR:    PlottingFiles/ContoursuftBimodal
%           Plotting file for uft bimodal experiments. Requires that the
%           values of Acel, uinterf, and temp are already computed, as well 
%           as uft. The idea of this example is to turn uft on for a frac
%           specific fraction of a day and then off for
%           the remaining fraction of the day. This will happen a set
%           amount of times, which we control. The average uft however
%           will always be 1. This is done by choosing the day fraction
%           such that the average is 1 when the 'on' value is 24.
%           Basically, to simplify it in a way, we choose dayfrac in such a
%           a way that the average over one hour is 24, which allows the
%           rest of the 23 hours to be zero in a day. 
%           There are optional plotting schemes for Pt and uft, if required. 
%           Lambda is plotted in all cases. 
%           The user is prompted for whether or not she would like to save
%           the plots. This option will then be applied to all following
%           plots. The user is then prompted to select which plot she would
%           like to see (and, if applicable, save). After it is ploted, the
%           user is asked if she would like to run the code again, in order
%           to plot or re-plot other figures. An input of 0 exits the code,
%           and finally lambda is plotted (along with other solely time
%           dependent variables such as P0t and uft).
% INPUT: 
%           No input variables
%          
% OUTPUT:   Main outcomes: 
%           Contours of one or some of:
%           Area, Temperature and Velocitu 
%           Lambda, in all cases. Where applicable, P0t and uft. 
% ADDITIONAL COMMENTS: 
%
% ASSOCIATED FUNCTIONS:
%           TimeDependentMOL, coupledPde, TimeDependenuft, coupledPdeuft,
%           and associated functions. 

%close all 
xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

% Here we can choose how many days we can display. In order to provide the
% best accuracy, we 
n = floor(T/d); % total number of days the simulation was run for

prompt = ' How many days do you want to display on the plot ? ';
k = input(prompt); 
step =(n-k)*d; % number of days we display on the plot is given by n*d - step
tnew = t - step;
indx = find(tnew>0,1);
indx = indx-1;
tvector = t*ones(1,2*K)-step; 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
tvectorint = t*ones(1,2*K+1)-step;
xvectorint =  [lam*xint',lam + (L-lam)*xint(2:end)'];

prompt = 'Do you wish to save data? (yes == 1) ';
sav = input(prompt);
r = 1; 
while r == 1
    prompt = [' Which plots do you want to produce? \n ', ...
        'Temperature == 1 \n ', 'Area == 2 \n ', 'Velocity == 3 \n ', ...
        'Note: in all cases we produce lambda, and where applicable, u(t) and P(t) \n '];
    plt = input(prompt);

    switch plt
        case 1 
            figure;
            title('Temperature')

            contourf(tvector(indx:end,:), xvector(indx:end,:),temp(indx:end,:),20,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            xlim([0 tvec(end)-step])
            maxtemp = max(max(temp))
            if sav==1
                axis off
                print(gcf, '-dpng', '-r600', '-painters', 'Temperature.png')
            end
        case 2 
            figure; 
            contourf(tvector(indx:end,:), xvector(indx:end,:),Acel(indx:end,:), 20,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            xlim([0 tvec(end)-step])
            if sav==1
                axis off
                print(gcf, '-dpng', '-r600', '-painters', 'A.png')
            end
        case 3
            
            %PLOTTING U
            figure; 
            contourf(tvectorint(indx:end,:), xvectorint(indx:end,:),uinterf(indx:end,:),80,'LineColor', 'none')
            xlim([0 tvec(end)-step])
            maxu = max(max(uinterf))
            caxis([min(min(uinterf)) 5])
            ax = gca;
            ax.YDir = 'reverse';
            if sav==1
                axis off
                print(gcf, '-dpng', '-r600', '-painters', 'Velocity.png')

            end
    end
    prompt = 'Do you want to run the code again? (yes == 1 ) \n ';
    r = input(prompt);
end
% PLOTTING lambda
figure; 
plot(t-step, lam);
xlim([0 t(end)-step])
title('lambda ')
xlabel('t')
if sav == 1 
   csvwrite('lam.csv',[t-step, lam]);
end
if uftval ==1
    figure;
    plot(t-step,uftvec);
    xlim([0 t(end)-step])
    title('uf')
    xlabel('t')
    if sav==1
        csvwrite('uft.csv',[t-step,uftvec]);
    end
end

if P0tval == 1
    figure
    plot(t-step,P)
    xlim([0 t(end)-step])
    title('P')
    xlabel('t')
    csvwrite('P.csv',[t-step, P]);
end
return

%code to calculate average P0


%indx1 = find(t==first(step),1);

indx1 = find(t>=step,1);

indx2 = N;

[tunique, index] = unique(tvec);

Pinterp = interp1(tunique,P(index),linspace(t(indx1),t(indx2),100));
mean(Pinterp)
Paverage = sum(P(indx1:indx2))/(indx2-indx1)

