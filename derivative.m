function d = derivative(f, dx)

d(1) = (f(2) - f(1)) / dx; % forward differences at node=1
d(length(f)) = ( f(end) - f(end-1) ) / dx; % backward differences at node=end

ndf = 2:(length(f)-1); 
d(ndf) = (f( ndf+1) - f(ndf-1)) / (2 * dx); % central differences for the middle

end