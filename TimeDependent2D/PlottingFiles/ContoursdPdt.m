% DATE:     2020 
% DESCR:    PlottingFiles/ContoursdPdt
%           Plotting file for Control systems. Requires that the
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



patt     = 5; 
prompt   = 'How many repetitions of the pattern length would you like to see ? \n';
k        = input(prompt);
rep      = length(first)/patt-k; 
step     = round(patt*rep);


% For pattern, the above step now starts at where the large addition is. In
% order to get t=0 to be somewhere 'reasonable', I need to shift another
% extra bit, which I choose to be the large t spacing / 3, i.e.

extrabit = 3*(first(step+1) - first(step))/4; 

%For the constant case,we want t0 to start at 0.0374, which is when the pattern
% case starts. So, we must add a value calculated specifically for the constant case.
% For pattern, remove this constant. For
%a different comparison, you must calculate the shift yourself.
shift = first(step) + 3*(first(step+1) - first(step)) + 0.0761;%+extrabit ;%- 0.0938; % The 0.0938 value is a fitting that I do to make
% the peak of the addition for 'day' 1, to match the one for CS2 . This is
% for figure 3.2.6 in my thesis, probably, hopefully
if m == 1
    
else
    pat_t0 = 0.0374;
    val = first-shift; 
    extrashift = pat_t0-val(find(val>0,1)); 
    shift = shift-extrashift; %+ add_shift; %-0.0194;
end


cylvec   = [first-shift, cyladd] ;
% c1 is the non-dimensional length of the cylinders. 
% here, cylvec is the number of cylinders added at each timestep
cylvec(:,2) = cylvec(:,2)/c1;
avg         = sum(cyladd)/length(first);
disp(['average is ', num2str(avg)]);

if sav == 1
    csvwrite('DStCyl.csv', cylvec); 
end


tnew = t - shift; 
indx = find(tnew>0,1); 
indx = indx -1;

xint       = linspace(0,1,K+1)';
xcel       = linspace(xint(2)/2,1-xint(2)/2,K)';
tvector    = t*ones(1,2*K)-shift; 
xvector    = [lam*xcel',lam + (L-lam)*xcel'];
tvectorint = t*ones(1,2*K+1)-shift;
xvectorint =  [lam*xint',lam + (L-lam)*xint(2:end)'];

r = 1; 
% PLOTTING THETA
while r == 1
    prompt = [' Which plots do you want to produce? \n ', ...
        'Temperature == 1 \n ', 'Area == 2 \n ', 'Velocity == 3 \n ', ...
        'Note: in all cases we produce lambda, and where applicable, u(t) and P(t) \n '];

    plt = input(prompt);

    switch plt
        case 1 
            figure;          
            contourf(tvector(indx:end,:), xvector(indx:end,:),temp(indx:end,:),20,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            xlim([0 tvec(end)-shift])
            if sav==1
                axis off
                print(gcf, '-dpng', '-r600', '-painters', 'Temperature.png')
            end

        case 2
            % PLOTTING A
            figure; 
            contourf(tvector(indx:end,:), xvector(indx:end,:),Acel(indx:end,:), 20,'LineColor', 'none')
            ax = gca;
            ax.YDir = 'reverse';
            xlim([0 tvec(end)-shift])
            if sav==1
                axis off
                print(gcf, '-dpng', '-r600', '-painters', 'A.png')
            end
        case 3
            % PLOTTING U
            figure; 
            contourf(tvectorint(indx:end,:), xvectorint(indx:end,:),uinterf(indx:end,:),20,'LineColor', 'none')
            xlim([0 tvec(end)-shift])
            if max(max(uinterf))>5
                maxu = 5
            else
                maxu = max(max(uinterf))
            end
            caxis([min(min(uinterf)) maxu])
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
plot(t-shift, lam);
xlim([0 t(end)-shift])
if sav ==1
     csvwrite('lam.csv',[t-shift, lam]);
end

hold on
if P0tval==0
    plot(t-shift,P);
     xlim([0 t(end)-shift])
    if sav==1
        csvwrite('P.csv',[t-shift, P])
    end
elseif P0tval==1
    plot(t,P0t(t));
    if sav==1
        csvwrite('P0t.csv',[t, P0t(t)])
    end
end
title('lambda and P')
xlabel('t')

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

