function val = uft2(t)
    global tend d 
    val = 0; 
    tstart = tend-d; 
    period = d/32; % 45 minutes 
    frac = (5/2700).*period; % 5 seconds of the day
    
    tsect = rem(t-tstart,period); 
    
    if tsect <= frac
        val = 1;
    end

end