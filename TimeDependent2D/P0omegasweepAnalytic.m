

% Find A0, \bar(u0), th0, and lam0 
% Solve the steady state 
close all

% Fix P0experiments to solve for sinusoidal Pin first
ParametersDefinition
%P0experiments % This gives Acel, temp, and uint, which have N rows and K columns
global N K Gamma P0 P0t    
% These are column vectors

% Change filename to matth what we want to import 
data = csvread('SSDataP01.csv');
A0 = data(2*K+1:3*K)'; 
u0bar = data(8*K+2:10*K+2)';
th0 =  data(4*K+1:6*K)';
    
lam0 = data(10*K+3); 
% We define the sinusoidal Pin, with respect to the constant P0 used
% for the steady state
Deltavalues = linspace(0,1,50)./P0;
lam1tilde = zeros(N,length(Deltavalues));
for i = 1:length(Deltavalues)
%% FOR SINUSOIDAL P
DeltaP = Deltavalues(i); 
n = 10;

P0t = @(t) P0 + DeltaP.*sin(t); % base case 

tau = linspace(0,2*pi,N)'; 
Pin = P0t(tau);
Pinbar = mean(P0t(tau));
Pintilde = Pin- Pinbar;
Pintau = Pintilde.*tau;
%DeltaP = 0.5;  % not sure about this
%omegaP = 5;
% top , from 0 to lam0
dx = 1/(K-1);
xvec = linspace(0,1,K);
xtop = lam0.*xvec; 
%bottom, from lam0 to 1

%fac2 = derivative(A0.*Cint,dx); 
%tau = omegaP.*t;

phitau2 = zeros(N,1);
phitauAn = phitau2;

%tau = linspace(0,2*pi,N)';
dt = 1/(N-1);
% analytic expression for Pintau

%vector expression for Pintau
term1 = (1/(2*pi))*trapz(Pintau).*((2*pi-0)*dt); 

phitau= term1 + cumtrapz(Pintilde)*(2*pi-0).*dt; 

%lam1tilde1 = -phitau*(1./(3*A0prime(K)*exp(-Gamma*th0(K))));
lam1tilde(:,i) = -phitau*(1./(3*A0primelam0*exp(-Gamma*th0(K))));

end

% % PLOTTING lam1tilde num for different delta values
 figure; 
% %tvector = t*ones(1,2*K); 
% % xvector = [lam0*xcel',lam0 + (L-lam0)*xcel'];
 taumat = tau*ones(1,length(Deltavalues));
% 
Deltamat = ones(N,1)*Deltavalues;
contourf(taumat, Deltamat, lam1tilde, 30,'LineColor', 'none')
colormap(jet)
ax = gca;
%ax.YDir = 'reverse';
%title('$\tilde{A}_1$ Analytic, $\omega_1=10$','Interpreter','latex','FontSize',14);
colorbar
%caxis([-0.4 0.2])
% if sav==1
%     axis off
%     colorbar off
%     print(gcf, '-dpng', '-r600', '-painters', 'SavedPlots/A1tildeAnalytic.png')
%    
% end
% 
