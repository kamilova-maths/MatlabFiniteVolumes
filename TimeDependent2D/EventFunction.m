function [value, isterminal, direction] = EventFunction(t,y)
% DATE:  2020
%
% DESCR:    [value, isterminal, direction] = EventFunction(t,y)
%           Event function for control system 1. This allows ode15s.m to
%           stop when P is below a certain tolerance (which I have set to
%           0.01) automatically. This is when the code itself decided when
%           to add the cylinders as opposed to us imposing it. 
%
% INPUT: 
%           t:  Particular t value using by the internal MATLAB routine
%           ode15s.m
%           y:  Long vector that includes all the variables of interest. It
%           is of the same length as y0, the initial condition vector,
%           provided in the main code that invokes ode15s.m
% OUTPUT:   value: this is an indicator value, which when 0, activates 
%           the event. Otherwise, it just carries on with ode15s.m as
%           usual. We set value = 0 when P is smaller than 0.01, the
%           tolerance we have decided to be tol_P.
%           isterminal: indicator value, which is set to 1 so that the
%           code is interrupted when event is activated
%           direction: I really don't know what this is. The official
%           description is: if all zeros are to be located (the default). 
%           A value of +1 locates only zeros where the event function is 
%           increasing, and -1 locates only zeros where the event function
%           is decreasing. Specify direction = [] to use the default value 
%           of 0 for all events. I set it to zero ... 
% ADDITIONAL COMMENTS: 
%           This code is called by coupledPdeWithdPdt.  
%
% ASSOCIATED FUNCTIONS:
%           ParametersDefinition : This is where all the parameters are
%           set, according to the specific need of the example.
%           bvpinit, bvp4c: to solve the ODEs
%           usolution: Routine that solves the u problem, which is
%           independent of time, so it does not need to be in the time
%           vector. We compute the form of u at each timestep t to use it
%           in the temperature and area equations. 
%           coupledPdeWithdPdt
global K

    if y(3*K+1)<0.01 % if P is less than this tolerance value, then indicator is zero.
        value  = 0; 
    else
        value  = 1; % otherwise, indicator is 1 
    end
    
    isterminal = 1; % Indicator value for stopping the function upon occurrence of the event 
    direction  = 0; %now what in the world could this be ? I truly have no idea. 
end