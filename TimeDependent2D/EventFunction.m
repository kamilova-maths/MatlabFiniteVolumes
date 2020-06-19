function [value, isterminal, direction] = EventFunction(t,y)
    global d dCdt
    value = mod(t,d);
    if value ==0
        dCdt = dCdt +1; 
    end
    isterminal = 0;
    direction = 0; 
end