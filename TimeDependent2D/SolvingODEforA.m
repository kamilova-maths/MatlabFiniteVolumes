function At= SolvingODEforA(t,A,n,mu,u,P)


dx = 1/(n-1);
for i=1:n
    if i==n
        At = P(i)/(3*mu) - u(i)*(0-A(i))/dx;
    else
        At = P(i)/(3*mu) - u(i)*((A(i+1)-A(i))/dx); 
    end
end


end