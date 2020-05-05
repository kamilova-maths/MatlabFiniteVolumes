%function [ th, A, u, x, t ] = TimeDependentFDfullMOL( th0, A0, D, gamma, P0, Pe, St, Bi, tha, T, L, K, N,uf,plt)
function yt = coupledPde(t,y)
% One vector to two vectors
global Pe Bi tha L K gamma P0 St D uf 
% for i=1:n
%     A(i) = y(i); 
%     th(i) = y(i+n);
% end
dx = L/K;
x = ( dx:dx:L )'; % Since we have Dirichlet boundary conditions, we don't need x=0

Atop = y(1:K); % This always imports the initial condition, so I don't really understand where to change u...
utop = y(1+K:2*K);
thtop = y(2*K+1:3*K); 

lam = y(3*K+1:4*K);

%Abot = y(4*K+1:5*K);
%ubot = y(5*K+1:6*K); 
thbot = y(4*K+1:5*K); 

% Both of these have to be the same value ... 
thmatch = thtop(end);
%thmatch = thbot(end); 

% viscosity
%mu = @(th) exp(-gamma*th); % Perhaps discretise this as well, but put a pin in that

% heat source
x1 = 5/7;
x2 = 6.5/7;
Qvalue = 1;
Q = Qvalue*(x>x1).*(x<x2);    % heat source  
%Qbot = Qvalue*(x>((1-x1)./lam)).*(x<((1-x2)./lam));    % heat source at the bottom

%% We solve everything for the TOP 
% We extract the necessary vectors 
A = Atop;
u = utop;
th = thtop;
%q = zeros(N,1);

% solve for u
temp1 = [0;0;th];
temp = (temp1(1:end-1)+temp1(2:end) ) /2; 
%tiph= 3*mu(temp); % evaluation 
tiph = 3*exp(-gamma*temp); 
tmpA = [ D; A  ];           % add ghost node to A

Dx2u = spdiags( [ tmpA(2:end).*tiph(2:end), -(tmpA(1:end-1).*tiph(1:end-1)+tmpA(2:end).*tiph(2:end)), tmpA(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], K, K );
Dx2u(1,2) = Dx2u(1,2) + tmpA(1).*tiph(1) / dx^2;              % include effect from Neumann BC
fu   = - St*(lam.^2).*(( tmpA(1:end-1) + tmpA(2:end) )/ 2);
fu(1) = fu(1) -2*tiph(1)*lam(1)*P0(1)/(3*dx); 
fu(end)= fu(end) - uf   * tmpA(end)* tiph(end)/dx^2;
u = Dx2u\fu;


tmpA =  (A + [A(2:end);1])/2; % extract A at the edges 
lamt = 1+(u(end)-u(end-1))./(tmpA(end)-tmpA(end-1));  % compute lamt with the edges


%tmpU = [ u; uf ];
tmpU = [ u; uf ] -lamt.*[x; 1] ; % size K+1 x 1 

% Source term

% add ghost node to A
tmpA = [D; A; 1]; % extend A by 1

%S  = -(A + [A(2:end);1]).*lamt./(2*lam); 
S = -A.*lamt./(lam); 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;
F  = ( F(1:end-1) - F(2:end) ) / dx; % F is the rhs to dA/dt
Arhs = F./lam+S; 

%% Solve for th at next time step (semi-implicit, parabolic)
%  % Assemble matrices
A = [Atop; ones(size(Atop))];
u = [utop; ones(size(utop))];
th = [thtop; thbot]; 

 % Assemble matrices
 tmpU = [ u; uf ]; 
% laplacian for theta
Dx2 = spdiags( ones(2*K,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], 2*K, 2*K );
Dx2(2*K,2*K-1) = 2/dx^2;

 % - first order derivative
temp = [D; A];         

g    = (2/Pe)./( A + [A(2:end);1] ) .* ( temp(2:end) - temp(1:end-1) )/dx - tmpU(2:end);  % notice u is treated explicitly
Dx1  = spdiags( (ones(2*K,1).*g).*[ 1, -1 ]/(2*dx), [-1,1], 2*K, 2*K )';
Dx1(2*K,2*K-1) = 0;
% - reaction term
Z = 2*Bi/Pe * spdiags( 1./sqrt(([A(2:end);1]+A)/2), 0, 2*K, 2*K );
    % - assemble matrix
M = ( Dx1 + 1/Pe * Dx2 - Z );
    % - assemble rhs
    xnew = linspace(0,1,2*K);
    Q = Qvalue.*(xnew>x1).*(xnew<x2);    % heat source  
fth = (2*(Bi/Pe) ./sqrt(([A(2:end);1]+A)/2)).*tha + Q';  

tht = M*th +fth;  % Then we impose the Dirichlet condition at X=1, which should match the

thttop = tht(1:K); thtbot = tht(K+1:end); 
% bottom 

% - assemble rhs

At = Arhs; 
ut = zeros(size(At));
% - solve
yt(1:K)= At;
yt(K+1:2*K) = ut ; 
yt(2*K+1:3*K) = thttop; 
yt(3*K+1:4*K) = lamt; 


% %% Then we solve everything for the BOTTOM
% % Abot = ones(size(A)); 
% % ubot = ones(size(A)); 
% 
% %tmpU = [ u; uf ];
% tmpU = lamt.*[x; 1]- 1 ; % size K+1 x 1 
% 
% g    = - (1./(1-lam)).*tmpU(2:end);  % notice u is treated explicitly
% Dx1  = spdiags( (ones(K,1).*g).*[ 1, -1 ]/(2*dx), [-1,1], K, K )';
% Dx1(2,1) = 0;
% % % - reaction term
% Z = 2*Bi/Pe * spdiags( ones(K,1), 0, K, K );
% %     % - assemble matrix
% M = ( Dx1 + 1/Pe * Dx2 - Z ); 
% %M(K,:) =0; M(:,K) = 0; 
%     % - assemble rhs
% fth = (Bi/Pe).*tha.*ones(K,1) + Qbot;  
% 
% fth(end) = fth(end) + thmatch.*(g(end)/(2*dx) + (1/(Pe*lam(end).^2.*(dx.^2)))); % we incorporate the matching condition 

 
%tht = M*thbot+fth;  
%At = zeros(size(tht)); 
%ut = zeros(size(At));
%ut = zeros(size(At)); 
% % - solve
% yt(4*K+1:5*K)= At;
% yt(5*K+1:6*K) = ut ; 
yt(4*K+1:5*K) = thtbot; 

yt = yt'; 

end

