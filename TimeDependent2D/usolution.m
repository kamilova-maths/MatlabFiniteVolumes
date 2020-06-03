function u = usolution(A,theta,lam,xmax)
% usolution  Computes solution to u, using A and theta. If the solution is
% computed separately from theta, then input a matrix of zeros. 
% INPUT:  A     : K x 1  vector for A
%         theta : K x 1 vector for theta - it could also be a constant of
%                zeros, for the cases where we do not calculate theta 
%         K     : discretisation in x 
% OUTPUT: u     : K x 1 vector for u 


global P0 K gamma uf St D
   
    %% Solve for u at next time step (laplacian)
%     temp = 3*mu([0;0;th(:,i)]);     % add ghost node to th
%     tiph = ( temp(1:end-1) + temp(2:end) ) / 2;
% u should only be computed here from 0 to 1 when it is in terms of X, and
% from 0 to L when it is in terms of x 
    dx = xmax/K;
    % Technically can do this in one line, but I can't figure out how -
    % this gives same results as above so we don't bother changing it
    %temp = [theta; theta(end)];
    %temp = (temp1(1:end-1)+temp1(2:end) ) /2; 
		theta = [0;theta];
    tiph = 3*exp(-gamma*theta);
    %tiph = [tiph; tiph(end)]; 
%     tiph = [tiph(1); tiph];
    % tiph = (tiph(1:end) + [tiph(2:end); tiph(end)])/2;
    
    Atmp  = [ D; A ];           % add ghost node to A
    % Add ghost node to A and ghost node to th (and extend by one term)
    Aint = ([ D; A] + [A;1] ) / 2;  

    Dx2u = spdiags( [ Atmp(2:end).*tiph(2:end), -(Atmp(1:end-1).*tiph(1:end-1)+Atmp(2:end).*tiph(2:end)), Atmp(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], K, K );
    Dx2u(1,2) = Dx2u(1,2) + Atmp(1).*tiph(1) / dx^2;              % include effect from Neumann BC
    %fu   = - St.*(lam.^2).*( Atmp(1:end-1) + Atmp(2:end) )/ 2;
    fu   = - St.*(lam.^2).*Aint(1:end-1);
    fu(1) = fu(1) - 2*(lam)*P0/(dx);  % include derivative (again, Neumann BC)
    
    fu(end)= fu(end) - uf   * Atmp(end)* tiph(end)/dx^2;
    u = Dx2u\fu;


end