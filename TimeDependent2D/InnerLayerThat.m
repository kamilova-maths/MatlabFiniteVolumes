
function yt = InnerLayerThat(t,y,X,dth0dx,dA0dx,u0hat0,th1hat0,ddth0dx,DeltaP)

global K Pe D Bi  tha Pinbar 
th = y(1:end);
dx = 1/(K-1);
%Pintildeatt = interp1(linspace(0,2*pi,N),Pintilde,rem(t,2*pi));
%thtmp = [bc; th; 0];% NOTE: I took the temp value from phi (continuity)
thtmp = [0; th; 0];% NOTE: I took the temp value from phi (continuity)
%Pintildeatt = DeltaP.*sin(t);
%A0tmp = [ D; A0 ];

Fw =-(1./(dx.*Pe)).*((thtmp(2:end)-thtmp(1:end-1))); 

Fw = ((Fw(1:end-1) - Fw(2:end) )./( dx)); 

% c1 = 2*Bi*tha/(Pe*sqrt(D)) + dth0dx(1)*dA0dx(1)/(Pe*D) - (dth0dx(1)/(Pe*D)).*Pinbar/(3*dA0dx(1));
% c2 = -u0hat0.*dth0dx(1)/(Pe*D);
% c3 = c1 + ddth0dx(1)/Pe;
% c4 = c2-th1hat0(1);
%S = c1 + c3.*Pintildeatt;
c3=-0;c4=0;
%S =  DeltaP.*exp(-X).*th1hat0(1).*(sin(t)+(1/Pe)*cos(t));
S =  DeltaP.*exp(-X).*th1hat0(1).*(sin(t)+(1/Pe)*cos(t));
yt = Fw + S;
%yt = yt'; 

end