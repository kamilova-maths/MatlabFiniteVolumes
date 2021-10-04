 

eta2 = @(tau,lam1) omega*Lam0 + lam1 - tau; 
tau_long = linspace(-100*pi, 100*pi, N)'; 

Fx3 = @(argtau,tau) dvex(0)*phiex(tau); % Fx for specific values of tau
    % first solve outside of the loop, so we can use the previous solution
    % as an initial guess at each point
tau_short = linspace(0, 2*pi, 200)';
lam1_val = zeros(size(tau_short));
Flam = @(lam1) duex(lam0).*lam1 - (1/3).*phiex(tau_short(1)) -  ...
           interp1(arg(tau_long), Fx(tau_long),eta2(tau_short(1), lam1),'spline');%-Fx3(omega*Lam0+lam1-tau_short(1),tau_short(1));
     % options = optimset('Display','off');
lam1_val(1) = fsolve(Flam, 0, options);
  
    for j = 2:length(tau_short)
        Flam = @(lam1) duex(lam0).*lam1 - (1/3).*phiex(tau_short(j))-...
              interp1(arg(tau_long), Fx(tau_long),eta2(tau_short(j), lam1),'spline');%-Fx3(omega*Lam0+lam1-tau_short(j),tau_short(j));
        lam1_val(j) = fsolve(Flam, lam1_val(j-1));
       % Flam(lam1_val(j));
    end
    
    %Flam = @(lam1) duex(lam0).*lam1 - (1/3).*phiex(tau_short)-Fx3(omega*Lam0+lam1-tau_short,tau_short);
    
figure;
plot(omega*t(indx:end)-2*pi*(n-1),lam1tildenum(indx:end))
hold on 
plot(tau_short, lam1_val)
xlim([0 2*pi])



% Just the correction by F
Fextra = interp1(arg(tau_long), Fx(tau_long),eta2(tau_short, lam1_val),'spline') ;

csvwrite('PlottingFiles/SavedPlots/lam1tildeNew.csv',[tau_short, lam1_val, Fextra]);

 