function d = derivative(f, dx)
% DATE:  2020
%
% DESCR:    d = derivative(f, dx)
%           Code that calculates the numerical derivative of a given function 
%	    f, given as a vector, with respect to a step size dx. Central 
%	    differences are used for the interior, with forward at the first node
%	    and backward at the last.  	
%
% INPUT: 
%           f:  vector of data points that describe some function f evaluated  at x
%           dx:  step size
% OUTPUT:   d: first derivative of f wrt x
% 
% ADDITIONAL COMMENTS: 
%      
%
% ASSOCIATED FUNCTIONS:
%           
%
%


d(1) = (f(2) - f(1)) / dx; % forward differences at node=1
d(length(f)) = ( f(end) - f(end-1) ) / dx; % backward differences at node=end

ndf = 2:(length(f)-1); 
d(ndf) = (f( ndf+1) - f(ndf-1)) / (2 * dx); % central differences for the middle

end
