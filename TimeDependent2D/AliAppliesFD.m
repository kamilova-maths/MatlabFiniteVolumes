function [ th ] = AliAppliesFD( th0, g, z, f, Pe, T, L, K, N )

dt = T/(K-1);
dx = L/N;

% Assemble matrices
% - laplacian
Dx2 = spdiags( ones(N,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], N, N );
Dx2(N,N-1) = 2/dx^2;
% - first order derivative
Dx1 = spdiags( (ones(N,1).*g).*[ 1, -1 ]/(2*dx), [-1,1], N, N )';
Dx1(N,N-1) = 0;


M = speye(N) - dt*( Dx1 + 1/Pe * Dx2 - speye( N ).*z' );

th = zeros(N,K);
th(:,1) = th0;

for i=1:K-1
    th(:,i+1) = M\( th(:,i) + dt*f(:,i+1) );
end


% % for plotting
% th = [ zeros(1,K); ...
%        th          ];
% 
% surf(th);

% % for debugging:
% th = t.*x.*(x-2);
% f  = x.*(x-2) - 2*g*t.*(x-1) -2*t/Pe + t.*x.*(x-2);

end

