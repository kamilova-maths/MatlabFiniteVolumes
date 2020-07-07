function u = usolution(A,theta,lam,xmax,P)
% usolution  Computes solution to u, using A and theta. If the solution is
% computed separately from theta, then input a matrix of zeros. 
% INPUT:  A     : K x 1  vector for A
%         theta : K x 1 vector for theta - it could also be a constant of
%                zeros, for the cases where we do not calculate theta 
%         K     : discretisation in x 
% OUTPUT: u     : K x 1 vector for u 


global K Gamma St D uf P0
    % I feel like I am missing one value of u right before this. 
    dx = xmax/(K);
	theta = [0;theta];
    tiph = 3*exp(-Gamma*theta);
    %Atmp = [D; A]; 
    Atmp = [2*D-A(1); A];
  
    % Add ghost node to A and ghost node to th (and extend by one term)
    Aint = ([D; A] + [A; 1]) /2; 
    
    % The order of this matrix is LOWER, DIAG, UPPER
    Dx2u = spdiags( [ Atmp(2:end).*tiph(2:end), -(Atmp(1:end-1).*tiph(1:end-1)+ ...
           Atmp(2:end).*tiph(2:end)), Atmp(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], K, K );
    %Dx2u(1,1) = Dx2u(1,1) + Atmp(1).*tiph(1) / dx^2; 
    Dx2u(1,2) = Dx2u(1,2) + Atmp(1).*tiph(1) / dx^2;   % include effect from Neumann BC

    fu   = - St.*(lam.^2).*Aint(1:end-1);
    fu(1) = fu(1) - 2*P*lam/(dx);  % include derivative (again, Neumann BC)
    %fu(1) = fu(1) - P*lam/(dx);  % include derivative (again, Neumann BC)
    fu(end)= fu(end) - uf* Atmp(end)* tiph(end)/dx^2;
    %fu(end)= fu(end) - (1./Aint(end))* Atmp(end-1)* tiph(end-1)/dx^2;
    u = Dx2u\fu;
end