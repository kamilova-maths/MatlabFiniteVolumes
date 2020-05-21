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
    temp = theta;
    %temp = [0;theta];
    %temp = (temp1(1:end-1)+temp1(2:end) ) /2; 
    tiph= 3*exp(-gamma*temp); 
    Atmp = A; 
    %Atmp  = [ D; A ];           % add ghost node to A
    % Atmp(end) = 1; 
    Dx2u = spdiags( [ Atmp(2:end).*tiph(2:end), -(Atmp(1:end-1).*tiph(1:end-1)+Atmp(2:end).*tiph(2:end)), Atmp(1:end-1).*tiph(1:end-1) ] / dx^2, [-1,0,1], K-1, K-1 );
    Dx2u(1,2) = Dx2u(1,2) + Atmp(1).*tiph(1) / dx^2;              % include effect from Neumann BC
    fu   = - St.*(lam.^2).*( Atmp(1:end-1) + Atmp(2:end) )/ 2;

    fu(1) = fu(1) -2*P0/(dx*lam(1));  % include derivative (again, Neumann BC)
    
    fu(end)= fu(end) - uf   * Atmp(end)* tiph(end)/dx^2;
    u = Dx2u\fu;
    
    u = [u; 1];


end