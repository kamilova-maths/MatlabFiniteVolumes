% DATE:     2020 
% DESCR:    PlottingFiles/ContoursAlternative
%           Plotting file for Control system 2. Requires that the
%           values of Acel, uinterf, and temp are already computed. There are
%           optional plotting schemes for Pt and uft, if required. Lambda
%           is plotted in all cases. 
%           The user is prompted for whether or not she would like to save
%           the plots. This option will then be applied to all following
%           plots. The user is then prompted to select which plot she would
%           like to see (and, if applicable, save). After it is plotted, the
%           user is asked if she would like to run the code again, in order
%           to plot or re-plot other figures. An input of 0 exits the code,
%           and finally lambda is plotted (along with other solely time
%           dependent variables such as P0t and uft).
%           Each pattern is typically of length 5, unless otherwise
%           specified. The values are chosen such that a comparison between
%           the constant control system and the pattern control system is
%           possible. If you require a different plotting, one that is
%           separate from this comparison, modify the code accordingly. In
%           particular, we must remove the constant value added after
%           'extrabit', to make the shift look nicer.
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

prompt   = 'Do you wish to save data? (yes == 1) ';
sav      = input(prompt);

prompt   = 'How many repetitions of the pattern length would you like to see ? \n';
rep        = input(prompt);


patt = 7; 
%rep  = 1; 
step = patt*rep*d;

tnew = t - step; 
indx = find(tnew>0,1); 
indx = indx -1;

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
%indx = N;
tvector = t*ones(1,2*K)-step; 
xvector = [lam*xcel',lam + (L-lam)*xcel'];
tvectorint = t*ones(1,2*K+1)-step;
xvectorint =  [lam*xint',lam + (L-lam)*xint(2:end)'];



%step = T-7*d;
%step = d;
if sav == 1
    list = patt*rep:floor(T/d);
    dvalues = d*list; 
    dvalues = dvalues-patt*rep*d; 
    cylvec = [dvalues(2:end)', cyladd(patt*rep:end)] ; 
    %cylvec = [dvalues(1:end)', cyladd(patt*rep:end)] ; 
    csvwrite('DStCyl.csv', cylvec); 
    njvalues = [dvalues(2:end)',cylvec(:,2)./onecyl];
    %njvalues = [dvalues(1:end)',cylvec(:,2)./onecyl];
    csvwrite('njvalues.csv',njvalues);
    csvwrite('height.csv',[njvalues(:,1), c1dim*njvalues(:,2)])
end

r = 1; 
% PLOTTING THETA
while r == 1
    prompt = [' Which plots do you want to produce? \n ', ...
        'Temperature == 1 \n ', 'Area == 2 \n ', 'Velocity == 3 \n ', ...
        'Note: in all cases we produce lambda, and where applicable, u(t) and P(t) \n '];

    plt = input(prompt);

    switch plt
        case 1 

            %PLOTTING THETA
            figure;

            contourf(tvector(indx:end,:), xvector(indx:end,:),temp(indx:end,:),20,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            xlim([0 tvec(end)-step])
            if sav==1
                axis off
                print(gcf, '-dpng', '-r600', '-painters', 'Temperature.png')
            end

        case 2 
            %PLOTTING A
            figure; 
            %tvector = t*ones(1,2*K)-step; 
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
            if max(max(uinterf))>5
                maxu = 5;
            else
                maxu = max(max(uinterf));
            end
            disp(['max u is ', num2str(maxu)])
            caxis([min(min(uinterf(indx:end,:))) maxu])
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
plot(t(indx:end)-step, lam(indx:end));
xlim([0 t(end)-step])

hold on
if P0tval==0
    plot(t(indx:end)-step,P(indx:end));
     xlim([0 t(end)-step])
    if sav==1
        csvwrite('P.csv',[t-step, P])
        csvwrite('lam.csv',[t-step, lam]);
    end
end
title('lambda and P')
xlabel('t')

if uftval ==1
    figure;
    plot(t(indx:end)-step,uftvec(indx:end));
    xlim([0 t(end)-step])
    title('uf')
    xlabel('t')
    if sav==1
        csvwrite('uft.csv',[tvecuft,uftvec]);
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

