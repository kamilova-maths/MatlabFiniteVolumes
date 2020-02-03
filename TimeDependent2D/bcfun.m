function res = bcfun(ya,yb,P0,R0,Tin)
res = [ ya(1,:)-P0;
        ya(2,:)-R0^2; 
        yb(3,:);
        ya(4,:)-Tin 
        ];
end