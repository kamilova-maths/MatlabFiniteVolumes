% DATE:     2020 
% DESCR:    LargeOmega1AsymptoticsIsothermal
%           Main code looking at the fluctuations in the isothermal case (Gamma = 0).
%           Uses TimeDependentMOLIsothermal, ParametersDefinition, 
%           coupledPdeIsothermal, and usolution to solve pdes. 
%           First we calculate an analytical solution, where we choose the
%           Stokes number such that lambda  = 1. Then we calculate u for
%           this Stokes number and compare and print the difference. The
%           discrepancy is attributed to discretisation errors, and there
%           is no way to circumvent it. Then we use either simple or steady
%           conditions to start the code. We impose P0t as a function, but
%           can be chosen to be a constant value if required. The results
%           are plotted at discrete timesteps to show convergence towards
%           steady state. 
% INPUT: 
%           No input variables
%          
% OUTPUT:   Main outcomes: 
%           Acel: The N x 2*K matrix storing the cell values for A 
%           uinterf: The N x (2*K +1) matrix storing all interface values for
%           u
%           
% ADDITIONAL COMMENTS: 
%
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           coupledPdeIsothermal: This is where the isothermal pdes are, that will be then
%           solved with method of lines with ode15s


close all 



redo = 1;
while redo == 1
    TimeDependentMOLIsothermal % This gives Acel, temp, and uinterf, which have N rows and K columns
    prompt = ' Would you like to refine this further? (yes == 1) \n [Remember to rename SSNEW] ';
    redo = input(prompt);
end

global N K L P0 P0t    

% Our steady state should match the analytical steady state obtained in 
% TimeDependentMOLIsothermal , and we have to transpose them so that they
% are row vectors
A0 = A0';
u0bar  = 1./A0; 
     

% We define the sinusoidal Pin, with respect to the constant P0 used
% for the steady state
%% FOR SINUSOIDAL P
xvec = linspace(0,1,K);
xtop = lam0.*xvec; 
tau      = linspace(0,2*pi,N)';
method = 'ana';
switch(method)
    case 'num'
        %tau      = linspace(0,2*pi,N)'; 
        
        Pin      = P0t(tau/omega);
        Pinbar   = mean(P0t(tau/omega));
        Pintilde = Pin- Pinbar;
        Pintau   = Pintilde.*tau;

        PinAn      = @(tau) P0+DeltaP*sin(tau); 
        PintauAn   = @(tau) DeltaP*sin(tau).*tau; 
        PintildeAn = @(tau) DeltaP*sin(tau);
        PinbarAn   = @(tau) P0; 
        % top , from 0 to lam0

     
      
        dx   = (xtop(end)-xtop(1))/(K-1);
        %bottom, from lam0 to 1
        %xbot = lam0 + (L-lam0).*xvec; % for the integrals, this is okay, as lam0 is a constant.
        % Integral from lam0 to x, where x goes from lam0 to 1 .The negative of
        % this is the integral we want, i.e. integral from x to lam0 where x goes
        % from lam0 to 1


        phitau2  = zeros(N,1);
        phitauAn = phitau2;

        %tau = linspace(0,2*pi,N)';
        dt = 1/(N-1);
        % analytic expression for Pintau
        term1An = (1/(2*pi))*integral(PintauAn,0,2*pi);

        %vector expression for Pintau
        term1   = (1/(2*pi))*trapz(Pintau).*((2*pi-0)*dt); 

        for i = 1:N
            phitauAn(i) = term1An + integral(PintildeAn,0,tau(i)); 
            phitau2(i)  = term1 + trapz(Pintilde(1:i))*(tau(i)-0)*(1/(i)); 
        end
        phitau = term1 + cumtrapz(Pintilde)*(2*pi-0).*dt; 
        check  = cumtrapz(Pintilde)*(2*pi-0).*dt; 

        fA0     = flip(A0(1:K)); 
        % The upper limit in this integral is 1 because dx already contains the
        % lam0 factor, remember that xtop goes from 0 to lam0, so dx is lam0-0/K,
        % and so the upper limit is 1*xtop(end)
        Cint    = flip(-cumtrapz(1./(fA0)).*(-1)*dx);  
        fac1    = smooth(derivative(A0(1:K).*Cint,dx));
        fac2    = smooth(derivative(A0(1:K),dx).*Cint + A0(1:K).*derivative(Cint,dx)); 
        A1tilde1 = -(1/3)*phitau*fac2';
        A1tilde = (1/3)*phitau*(1-derivative(A0(1:K),dx).*Cint);


        %A0prime = derivative(A0,dx); 
        A0primelam0 = (St.*trapz([A0(1:K), 1]).*(1).*dx + Pinbar)./(3*(1)); 

        %lam1tilde1 = -phitau*(1./(3*A0prime(K)*exp(-Gamma*th0(K))));
        lam1tilde = -phitau*(1./(3*A0primelam0));


        u0tilde = (Pintilde./3)*Cint(1:K);

       
    case 'ana'


%% Calculating the analytical solution of the 'extra bit' for Acal

   % x     = linspace(0,lam0,K);
    F   = sqrt(6*St*D-P0^2);
    a1 = atan(P0/F); 
    Pstar =  6*a1/(F); 
    Aex   = @(x)(F^2/(6*St)).*tan(F.*(x+Pstar)./6).^2 - P0.^2./(6*St)+D ; 
    uex = @(x) 1./Aex(x); 
    c0 = -DeltaP; 

    phiex = @(tau) -DeltaP*cos(tau) ; 
    % Now we define the derivatives explicitly to prevent approximation errors
    % due to finite differences

    dAex = @(x) F^3.*(sec((F*x)./6 + a1).^2).* tan((F.*x)/6 + a1)./(18*St); 

    ddAex = @(x)-((F^4 .*(-2 +cos((F.*x)./3 + 2.*a1)).*(sec((F.*x)/6 + a1).^4))./(108*St)); 


    psiex = @(x) (-P0 + F.*tan((F.*x)./6 + a1 ))./St ;

    uhatex = @(x) (St.*(F.*(-x + lam0) - 3*sin((F.*x)./3 + 2.*a1) +  ...
    3.*sin((F.*lam0)./3 + 2.*a1)))./F^3 ;

    duhatex = @(x)(St.*(-F - F.*cos((F.*x)./3 + 2 .*a1)))./F^3 ;
    c1 = (uex(0)/uhatex(0))*(dAex(0)*uhatex(0) - (1/3)); 
    
    Xhat =  @(x) omega.*psiex(x); 
    vex = @(x) uhatex(x)./uex(x); 
    dvex  = @(x) -(1/36).*(sec(a1 + (F.*x)./6).^3).*(12.*cos(a1 + (F.*x)./6) - ... 
   3.*cos(a1 - (1/6).*F.*(x - 2.*lam0)) + 3.*cos(3.*a1 + (1/6).* F.*(x + 2.*lam0)) + ... 
   2.*F.*(x -lam0).*sin(a1 + (F.*x)./6)); 

    A1tilde_fun = @(x,tau) -phiex(tau)*dvex(x); 
    A1tilde = A1tilde_fun(xtop,tau) ;
     
      % Extra plots for Acal
    %figure; 
    arg = @(tau) -tau - D*uhatex(0)*phiex(tau);
    %v0ex0 = -(D*St/(3*F^2))*(6-P0*lam0-(1/F)*(3*P0*sin((1/3)*F*lam0 + 2*atan(P0/F)))) ;

    Fx = @(tau) dvex(0)*phiex(tau);
    argtau = arg(tau);
    Fx2 = @(argtau) dvex(0)*phiex(tau);
    %plot(arg(tau),Fx(tau)); 
    tau2 = linspace(-2*pi,4*pi,300)';
    valuesmatrix=[arg(tau2),Fx(tau2)];
    csvwrite('PlottingFiles/SavedPlots/FetaDelta1.csv', valuesmatrix)
end
tval = t(end)-2*pi/omega; 
indx = find(t>=tval,1);
indx = indx -1; 

xint = linspace(0,1,K+1)';
xcel = linspace(xint(2)/2,1-xint(2)/2,K)';

xmatrix    = [lam*xcel',lam + (L-lam)*xcel'];
xmatrixint =  [lam*xint',lam + (L-lam)*xint(2:end)'];
xmat2      = lam*xcel';
xmat1int   = lam*xint';


Anum     = zeros(N,K);
u0barnum = zeros(N,K); 
for i = 1:N
    Anum(i,:) = interp1(xmatrix(i,:),Acel(i,:),linspace(0,L,K),'pchip'); 
end

xforu = linspace(0,1,K+1)*lam0;
for i = 1:N
     u0barnum(i,:) = interp1(xmatrixint(i,:),uinterf(i,:),linspace(0,L,K),'pchip'); 
end


omegat      = omega*t(indx:end);

AAvg        = mean(Anum(indx:end,:));

u0Avg       = mean(u0barnum(indx:end,:));


u0tildenum  = u0barnum - ones(N,1)*u0Avg; 


A1tildenum  = omega*(Anum-ones(N,1)*AAvg);

%A1tildenum2 = smooth(A1tildenum);

lam1tildenum = omega*(lam-mean(lam(indx:end))); 


run('PlottingFiles/ContoursLargeOmega1')

% Subtracting A1tilde from A1tildenum to isolate just the oscillations
omegat = omega*t - 2*pi*(n-1); 
A1tilde_small = A1tilde_fun(xtop,omegat(indx:end)); 
A1osc = A1tildenum(indx:end,:)-A1tilde_small;
%omegatmat_small = omegat(indx:end)*ones(1,K);  

% Option 3 - matrix form 
eta = @(x,tau) Xhat(x) - tau - vex(x).*phiex(tau);
tau_long = linspace(-100*pi, 100*pi, N)'; 
Fmat = interp1(arg(tau_long), Fx(tau_long),eta(xtop, omegat(indx:end)));

figure; 
contourf(omegatmat(indx:end,:), xmat(indx:end,:), Fmat, 100,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
xlim([0,2*pi])
ylim([0 lam(end)])
xlabel('$\tau$', 'Interpreter', 'latex')
ylabel('$x$', 'Interpreter', 'latex')
colorbar
set(gca, 'TickLabelInterpreter', 'latex');
if sav==1
    axis off
    colorbar off
    print(gcf, '-dpng', '-r300', '-painters', 'PlottingFiles/SavedPlots/FmatOmega100Dp5.png')

end


figure; 
contourf(omegatmat(indx:end,:), xmat(indx:end,:), A1osc, 100,'LineColor', 'none')
ax = gca;
ax.YDir = 'reverse';
xlim([0,2*pi])
ylim([0 lam(end)])
caxis([min(min(Fmat)), max(max(Fmat))])
if sav==1
    axis off
    colorbar off
    print(gcf, '-dpng', '-r300', '-painters', 'PlottingFiles/SavedPlots/A1OscOmega100Dp5.png')

end
%G1Avg = mean(A1osc);
%u1Avg = u0Avg-uex(x);
%q = Aex(x).*u1Avg + uex(x).*G1Avg;
%
% 
% Lam0 = psiex(lam0);
% fun = @(tau, lam1) -dAex(lam0)*lam1 +duhatex(lam0)*phiex(tau) - interp1(arg(tau_long), Fx(tau_long), omega*Lam0 + lam1 -tau); 
% lam1_tau = zeros(size(tau2));
% lam1_tau_prev = 0;
% tau2 = linspace(0, 2*pi, 50)';
% for i = 1:length(tau2)
%     tau_now = tau2(i);
%     fun_zero = @(lam1) fun(tau_now,lam1);
%     lam1_tau(i) = fsolve(fun_zero,lam1_tau_prev);
%     lam1_tau_prev = lam1_tau(i);
% end
% 
% 
