function u = usolution(A,theta,lam,xmax,P0t)
% usolution  Computes solution to u, using A and theta. If the solution is
% computed separately from theta, then input a matrix of zeros. 
% INPUT:  A     : K x 1  vector for A
%         theta : K x 1 vector for theta - it could also be a constant of
%                zeros, for the cases where we do not calculate theta 
%         K     : discretisation in x 
% OUTPUT: u     : K x 1 vector for u 


global K Gamma uf St D
   
 
    dx = xmax/K;

	theta = [0;theta];
    tiph = 3*exp(-Gamma*theta);

    
    Atmp  = [ D; A ];           % add ghost cell to A
    % Add ghost node to A and ghost node to th (and extend by one term)
    Aint = ([ D; A] + [A;1] ) / 2;  
    % The order of this matrix is LOWER, DIAG, UPPER
    Dx2u = spdiags( [ Atmp(2:end).*tiph(2:end), -(Atmp(1:end-1).*tiph(1:end-1)+ ...
        Atmp(2:end).*tiph(2:end)), Atmp(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], K, K );
    %Dx2u(1,1) = Dx2u(1,1) + Atmp(1).*tiph(1) / dx^2; 
    Dx2u(1,2) = Dx2u(1,2) + Atmp(1).*tiph(1) / dx^2;              % include effect from Neumann BC
    %fu   = - St.*(lam.^2).*( Atmp(1:end-1) + Atmp(2:end) )/ 2;
    %fu   = - St.*(lam.^2).*Aint(2:end);
    fu   = - St.*(lam.^2).*Aint(1:end-1);
    fu(1) = fu(1) - 2*P0t*lam/(dx);  % include derivative (again, Neumann BC)
    %fu(1) = fu(1) - P0t*lam/(dx);  % include derivative (again, Neumann BC)
    
    %fu(end) = fu(end) - uf * 3*exp(-gamma*phi0) / dx^2; 
    
    %fu(end)= fu(end) - uf   * Atmp(end-1)* tiph(end-1)/dx^2;
    fu(end)= fu(end) - uf   * Atmp(end)* tiph(end)/dx^2;
    u = Dx2u\fu;


end