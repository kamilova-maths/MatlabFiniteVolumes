function [value, isterminal, direction] = EventFunction(t,y)
global K

    if y(3*K+1)<0.01
        value = 0; 
    else
        value = 1;
    end
    isterminal = 1;
    direction = 0; %now what in the world could this be ?
end