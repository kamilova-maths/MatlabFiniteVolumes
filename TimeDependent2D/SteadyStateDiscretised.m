
% THIS CODE DOES NOT WORK - DO NOT USE THIS CODE - EVER EVER EVER


function [P, A, th] = SteadyStateDiscretised(N,gamma,Qvalue,x1,x2,eps,St,tha,Bi,Pe,P0,D,L,H)

% initialisation
th = zeros(N,1);
A  = zeros(N,1);
P = zeros(N,1); 

A(1,1) =D;
th(1,1)=0;
P(1,1)=P0;

% domain
dx = L/N;

x = ( dx:dx:L )';

% viscosity
mu = @(th) exp(-gamma*th);

% heat source
 Q = Qvalue*(x>x1).*(x<x2);    % heat source  
%q = zeros(N,1);


% laplacian for theta
Dx2 = spdiags( ones(N,1).*[ 1, -2, 1 ]/dx^2, [-1,0,1], N, N );
Dx2(N,N-1) = 2/dx^2;






for i=2:N
    %% Solve for P at next step (explicit)
    P(i)=St.*A(i-1)*dx+P(i-1);
    A(i)=((P(i-1).*A(i-1))./(3.*mu(th(i-1)))).*dx+A(i-1);
        
    
    if i==2
        
        th(i) = (1/(A(i)-dx.*Pe)).*(2.*sqrt(A(i-1)).*Bi.*dx^2*(-tha)-A(i-1).*(dx^2.*Pe.*Q(i-1)));    
    elseif i==N
        th(i) = th(i-1);
    else
        prefac= 1/(A(i)-dx.*Pe);
        thaterm= 2.*sqrt(A(i-1).*Bi.*dx^2).*(th(i-1)-tha);
        midterm= (A(i)-Pe.*dx).*th(i-1);
        anterm= A(i-1).*(dx^2.*Pe.*Q(i-1)-th(i-1)+th(i-2)); 
        th(i)= prefac.*(thaterm+midterm+anterm); 
    end
end  
 

end

