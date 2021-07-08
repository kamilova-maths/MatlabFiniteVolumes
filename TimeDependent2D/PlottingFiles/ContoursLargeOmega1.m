% DATE:     2020 
% DESCR:    PlottingFiles/ContoursLargeOmega1
%           File to plot contours for results where we impose a sinusoidal
%           P variation and we look at large frequency cases. Requires that the
%           values of Acel, uinterf, and temp are already computed. There are
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
taumat = tau*ones(1,K);
xmat1 = ones(N,1)*xtop;
omegatnew = omega*t-2*pi*(n-1);
indx = find(omegatnew>0,1);
indx = indx - 1; 
omegatmat = (omega*t-2*pi*(n-1))*ones(1,K);
omegat = omega*t - 2*pi*(n-1); 
xmat = ones(N,1)*linspace(0,L,K); 

prompt = 'Do you wish to save data? (yes == 1) ';
sav = input(prompt);
r = 1; 
while r == 1

prompt = [' Which plots do you want to produce? \n ', ...
        'Temperature == 1 \n ', 'Area == 2 \n ', 'Velocity == 3 \n ', ...
        'Everything else == 4 \n '];

    plt = input(prompt);
    prompt = ' Would you like to fix the axes ? (yes == 1 )';
    an = input(prompt);
    if an == 1
        prompt = ' With respect to analytical ( input == 0 ) or numerical (input == 1) result?' ;
        resc = input(prompt);
    end
    switch plt
        case 1 
            % ANALYTICAL
            figure; 
            contourf(taumat, xmat1, t1tilde(:,1:K), 100,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            xlim([0,2*pi])
            ylim([0 lam0])
            if an == 1
                if resc == 0
                    a = min(min(t1tilde));
                    b = max(max(t1tilde));
                else
                    a = min(min(t1tildenum(indx:end,:)));
                    b = max(max(t1tildenum(indx:end,:))); 
                end
            else
                a = min(min(t1tilde));
                b = max(max(t1tilde));
            end
            caxis([a, b])
            disp(['min theta is ', num2str(a) , '  ', 'max theta is ', num2str(b)])
            colorbar
            if sav==1
                axis off
                colorbar off
                print(gcf, '-dpng', '-r300', '-painters', 'SavedPlots/t1tildeAnalytic.png')

            end
            % NUMERICAL
            figure; 

            contourf(omegatmat(indx:end,:),xmat(indx:end,:),t1tildenum(indx:end,:),100,'LineColor','none')

            ax = gca;
            ax.YDir = 'reverse';
            xlim([0,2*pi])
            ylim([0, lam0])
            if an == 1
                if resc == 0
                    a = min(min(t1tilde));
                    b = max(max(t1tilde));
                else
                    a = min(min(t1tildenum(indx:end,:)));
                    b = max(max(t1tildenum(indx:end,:))); 
                end
            else
                 a = min(min(t1tildenum(indx:end,:)));
                 b = max(max(t1tildenum(indx:end,:))); 
            end
            caxis([a, b])
            disp(['min theta is ', num2str(a) , '  ', 'max theta is ', num2str(b)])
            colorbar
            if sav==1
                axis off
                colorbar off
                print(gcf, '-dpng', '-r300', '-painters', 'SavedPlots/t1tildeNumeric.png')

            end
        case 2
            % ANALYTICAL
            figure; 
            contourf(taumat, xmat1, A1tilde, 100,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            xlim([0,2*pi])
            ylim([0 lam0])
            if an == 1
                if resc == 0
                    a = min(min(A1tilde));
                    b = max(max(A1tilde));
                else
                    a = min(min(A1tildenum(indx:end,:)));
                    b = max(max(A1tildenum(indx:end,:))); 
                end
            else
                    a = min(min(A1tilde));
                    b = max(max(A1tilde));
            end
            caxis([a, b])
            disp(['min A is ', num2str(a) , '  ', 'max A is ', num2str(b)])
            colorbar
            if sav==1
                axis off
                colorbar off
                print(gcf, '-dpng', '-r300', '-painters', 'SavedPlots/A1tildeAnalytic.png')

            end

            % NUMERICAL 
            figure; 

            contourf(omegatmat(indx:end,:),xmat(indx:end,:),A1tildenum(indx:end,:),50,'LineColor','none')
            xlim([0,2*pi])
            ylim([0 lam0])
            ax = gca;
            ax.YDir = 'reverse';
            if an == 1
                if resc == 0
                    a = min(min(A1tilde));
                    b = max(max(A1tilde));
                else
                    a = min(min(A1tildenum(indx:end,:)));
                    b = max(max(A1tildenum(indx:end,:))); 
                end
            else
                    a = min(min(A1tildenum(indx:end,:)));
                    b = max(max(A1tildenum(indx:end,:))); 
            end 
            caxis([a, b])
            disp(['min A is ', num2str(a) , '  ', 'max A is ', num2str(b)])
            colorbar
            xlim([0,2*pi])
            ylim([0 lam0])
            if sav==1
                axis off
                colorbar off
                print(gcf, '-dpng', '-r300', '-painters', 'SavedPlots/A1tildeNumeric.png')

            end

        case 3

            % ANALYTICAL
            figure; 
            contourf(taumat, xmat1, u0tilde, 100,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            xlim([0,2*pi])
            ylim([0 lam0])
            if an == 1
                if resc == 0
                    a = min(min(u0tilde));
                    b = max(max(u0tilde));
                else
                    a = min(min(u0tildenum(indx:end,:)));
                    b = max(max(u0tildenum(indx:end,:))); 
                end
            else
                a = min(min(u0tilde));
                b = max(max(u0tilde));
            end 
            caxis([a, b])
            disp(['min u is ', num2str(a) , '  ', 'max u is ', num2str(b)])
            colorbar
            if sav==1
               axis off
               colorbar off
                print(gcf, '-dpng', '-r300', '-painters', 'SavedPlots/u0tildeAnalytic.png')

            end

            % NUMERICAL
            figure; 
            contourf(omegatmat(indx:end,:), xmat(indx:end,:), u0tildenum(indx:end,:), 100,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            xlim([0,2*pi])
            ylim([0 lam0])  
            if an == 1
                if resc == 0
                    a = min(min(u0tilde));
                    b = max(max(u0tilde));
                else
                    a = min(min(u0tildenum(indx:end,:)));
                    b = max(max(u0tildenum(indx:end,:))); 
                end
            else
                a = min(min(u0tildenum(indx:end,:)));
                b = max(max(u0tildenum(indx:end,:))); 
            end 
            caxis([a, b])
            disp(['min u is ', num2str(a) , '  ', 'max u is ', num2str(b)])
            colorbar
            if sav==1
                axis off
                colorbar off
                print(gcf, '-dpng', '-r300', '-painters', 'SavedPlots/u0tildeNumeric.png')

            end
        case 4

            lam1tildenum = omega*(lam-mean(lam(indx:end))); 
            lam1tildenum =lam1tildenum-(1/omega)*mean(lam1tildenum); 
            figure;
            plot(omega*t(indx:end)-2*pi*(n-1),lam1tildenum(indx:end))
            xlim([0 2*pi])
            hold on 
            plot(tau,lam1tilde)
            if sav == 1
                csvwrite('SavedPlots/lam1tilde.csv',[tau(1:16:end),lam1tilde(1:16:end)]);
                %csvwrite('SavedPlots/lam1tilde.csv',[tau,lam1tilde]);
                csvwrite('SavedPlots/phitau.csv',[tau,phitau]);
                csvwrite('SavedPlots/Pintilde.csv',[tau,Pintilde]);
                csvwrite('SavedPlots/lam1tildeNum.csv',[(omega*t(indx:end)-2*pi*(n-1)), lam1tildenum(indx:end)])
            end
    end
    prompt = 'Do you want to run the code again? (yes == 1 ) \n ';
    r = input(prompt);
end        
            

return
