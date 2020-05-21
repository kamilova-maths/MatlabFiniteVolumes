function res = bcfun(ya,yb,Tin)
global P0 D
res = [ ya(1,:)-P0;
        ya(2,:)-D; 
        yb(3,:);
        ya(4,:)-Tin 
        ];
end