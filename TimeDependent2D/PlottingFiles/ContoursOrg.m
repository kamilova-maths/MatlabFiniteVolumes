% DATE:     2020 
% DESCR:    PlottingFiles/ContoursOrg
%           Original file to plot contours for results. Requires that the
%           values of Acel, uint, and temp are already computed. There are
%           optional plotting schemes for Pt and uft, if required. Lambda
%           is plotted in all cases. 
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


% Defining all x values to be used
% x interfaces, from 0 to lamda0, rescaled so that it is 0 to 1
xint = linspace(0,1,K+1)';
% x cells, from 0 to lambda0, rescaled so that it is 0 to 1
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';
%indx = N;
tvector = t*ones(1,2*K); 
xvector = [lam*xcel',lam + (L-lam)*xcel'];

tvectorint = t*ones(1,2*K+1);
xvectorint =  [lam*xint',lam + (L-lam)*xint(2:end)'];

prompt = 'Do you wish to save data? (yes == 1) ';
sav = input(prompt);
r = 1; 
% PLOTTING THETA
while r == 1
    prompt = [' Which plots do you want to produce? \n ', ...
        'Temperature == 1 \n ', 'Area == 2 \n ', 'Velocity == 3 \n ', ...
        'Everything else == 4 \n ' ];

    plt = input(prompt);

    switch plt
        case 1 

            % % Save data to file
            figure;
            title('Temperature')

            xvector = [lam*xcel',lam + (L-lam)*xcel'];
            contourf(tvector, xvector,temp,20,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            if sav==1
                axis off
                print(gcf, '-dpng', '-r600', '-painters', 'Temperature.png')
            end
        case 2

            % PLOTTING A
            figure;
            contourf(tvector, xvector,Acel, 20,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            if sav==1
                axis off
                print(gcf, '-dpng', '-r600', '-painters', 'A.png')
            end

        case 3

            % PLOTTING U
            figure; 

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

            end
        case 4
            % PLOTTING lambda
            figure; 
            plot(t, lam);
            if sav == 1
               csvwrite('lam.csv',[t, lam]);
            end
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
    end
    
    prompt = 'Do you want to run the code again? (yes == 1 ) \n ';
    r = input(prompt);
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

