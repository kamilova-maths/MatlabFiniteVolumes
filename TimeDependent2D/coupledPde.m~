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

Abot = y(4*K+1:5*K);
ubot = y(5*K+1:6*K); 
thbot = y(6*K+1:7*K); 


% viscosity
%mu = @(th) exp(-gamma*th); % Perhaps discretise this as well, but put a pin in that

% heat source
x1 = 5/7;
x2 = 6.5/7;
Qvalue = 1;
Qtop = Qvalue*((x/lam)>x1).*((x/lam)<x2);    % heat source at the top
Qbot = Qvalue*(x>((1-x1)/lam)).*(x<((1-x2)/lam));    % heat source at the bottom

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




%   %% Solve for th at next time step (semi-implicit, parabolic)
%  % Assemble matrices
%  
% % laplacian for theta
Dx2 = spdiags( (1./lam.^2).*ones(K,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], K, K );
Dx2(K,K-1) = (lam(K).^2).*2/dx^2;
% 
%  % - first order derivative
 temp = [D; A];         
% 
g    = (2/Pe)./((lam.^2).*( A + [A(2:end);1] )) .* ( temp(2:end) - temp(1:end-1) )/dx - (1./lam).*tmpU(2:end);  % notice u is treated explicitly
Dx1  = spdiags( (ones(K,1).*g).*[ 1, -1 ]/(2*dx), [-1,1], K, K )';
%Dx1(K,K-1) = 0;
% % - reaction term
Z = 2*Bi/Pe * spdiags( 1./sqrt(([A(2:end);1]+A)/2), 0, K, K );
%     % - assemble matrix
M = ( Dx1 + 1/Pe * Dx2 - Z ); 

% - assemble rhs
fth = (2*(Bi/Pe) ./sqrt(([A(2:end);1]+A)/2)).*tha + Qtop;  

At = Arhs; 
tht = M*thtop +fth;  
ut = zeros(size(At));
% - solve
yt(1:K)= At;
yt(K+1:2*K) = ut ; 
yt(2*K+1:3*K) = tht; 
yt(3*K+1:4*K) = lamt; 


%% Then we solve everything for the BOTTOM
A = Abot;
u = ubot;
th = thbot; 

% solve for u
temp1 = [thtop(end);thtop(end);th];
temp = (temp1(1:end-1)+temp1(2:end) ) /2; 
%tiph= 3*mu(temp); % evaluation 
tiph = 3*exp(-gamma*temp); 
tmpA = [ 1; A  ];           % add ghost node to A


Dx2u = spdiags( [ tmpA(2:end).*tiph(2:end), -(tmpA(1:end-1).*tiph(1:end-1)+tmpA(2:end).*tiph(2:end)), tmpA(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], K, K );
Dx2u(1,:) =0;
%Dx2u(1,2) = Dx2u(1,2) + tmpA(1).*tiph(1) / dx^2;              % include effect from Neumann BC
fu   = - St*((1-lam).^2).*(( tmpA(1:end-1) + tmpA(2:end) )/ 2);
%fu(1) = fu(1) -1; 
%fu(end)= fu(end) - uf   * tmpA(end)* tiph(end)/dx^2;
u = Dx2u\fu;
u(1) = 1;  % Not sure at all that this is going to work xD


tmpA =  (A + [A(2:end);1])/2; % extract A at the edges 
lamt = 1+(u(end)-u(end-1))./(tmpA(end)-tmpA(end-1));  % compute lamt with the edges


%tmpU = [ u; uf ];
tmpU = lamt.*[x; 1]- [1; u] ; % size K+1 x 1 


% Source term

% add ghost node to A
tmpA = [D; A; 1]; % extend A by 1

%S  = -(A + [A(2:end);1]).*lamt./(2*lam); 
S = -A.*lamt./(1-lam); 
FL = tmpA(1:end-1).*tmpU;
FR = tmpA(2:end  ).*tmpU;
F  = ( FL+FR + abs(tmpU).*( tmpA(1:end-1) - tmpA(2:end) ) ) / 2 ;
F  = ( F(1:end-1) - F(2:end) ) / dx; % F is the rhs to dA/dt
Arhs = F./(1-lam)+S; 

%   %% Solve for th at next time step (semi-implicit, parabolic)
%  % Assemble matrices
%  
% % laplacian for theta
Dx2 = spdiags( (1./((1-lam).^2)).*ones(K,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], K, K );
Dx2(K,K-1) = ((1-lam(K)).^2).*2/dx^2;
% 
%  % - first order derivative
 temp = [D; A];         
% 
g    = (2/Pe)./(((1-lam).^2).*( A + [A(2:end);1] )) .* ( temp(2:end) - temp(1:end-1) )/dx - (1./(1-lam)).*tmpU(2:end);  % notice u is treated explicitly
Dx1  = spdiags( (ones(K,1).*g).*[ 1, -1 ]/(2*dx), [-1,1], K, K )';
Dx1(K,K-1) = 0;
% % - reaction term
Z = 2*Bi/Pe * spdiags( 1./sqrt(([A(2:end);1]+A)/2), 0, K, K );
%     % - assemble matrix
M = ( Dx1 + 1/Pe * Dx2 - Z ); 
    % - assemble rhs
fth = (2*(Bi/Pe) ./sqrt(([A(2:end);1]+A)/2)).*tha + Qbot;  

At = Arhs; 
tht = M*thtop +fth;  
ut = zeros(size(At));
% - solve
yt(1:K)= At;
yt(K+1:2*K) = ut ; 
yt(2*K+1:3*K) = tht; 
yt(3*K+1:4*K) = lamt; 









yt = yt'; 

end

