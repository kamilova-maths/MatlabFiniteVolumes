function [value, isterminal, direction] = EventFunctionIsothermal(t,y)
global K

    if y(K+1)<0.01
        value = 0; 
    else
        value = 1;
    end
    isterminal = 1;
    direction = 0; %now what in the world could this be ?
end