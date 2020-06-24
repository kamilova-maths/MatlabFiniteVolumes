function [value, isterminal, direction] = EventFunction(t,y)
    global Ldim uc first d
    %v1 = mod(t*Ldim/(86400*uc),1);
    v1 = mod(t,d); 
    tol = 1e-3;
    if abs(first-t)<tol
        value = 1; 
    elseif (abs(1-v1)<tol || abs(v1)<tol) && t>0.1
       value = 0;
        else
       value = 1;
    end
    isterminal = 1;
    direction = 0; %now what in the world could this be ?
end