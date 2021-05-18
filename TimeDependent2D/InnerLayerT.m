
function yt = InnerLayerT(t,y,dth0dx,dA0dx,u0hat0,th1hat0,ddth0dx)

global K Pe DeltaP
%D Bi tha Pinbar DeltaP
th = y(1:end);
dx = 1/(K-1);

%thend = interp1(linspace(0,2*pi,N),t1tilde,t);
%Pintildeatt = interp1(linspace(0,2*pi,N),Pintilde,rem(t,2*pi));
%Pintildeatt =  DeltaP.*sin(t);
phitauatt = DeltaP.*cos(t);
%phitauatt = interp1(linspace(0,2*pi,N),phitau,rem(t,2*pi));  
bc = th1hat0(1)*phitauatt; 
thtmp = [bc; th; 0];% NOTE: I took the temp value from phi (continuity)
%thtmp = [0; th; 0];% NOTE: I took the temp value from phi (continuity)

%A0tmp = [ D; A0 ];

Fw =-(1./(dx.*Pe)).*((thtmp(2:end)-thtmp(1:end-1))); 

Fw = ((Fw(1:end-1) - Fw(2:end) )./( dx)); 

% c1 = 2*Bi*tha/(Pe*sqrt(D)) + dth0dx(1)*dA0dx(1)/(Pe*D) - (dth0dx(1)/(Pe*D)).*Pinbar/(3*dA0dx(1));
% c2 = -u0hat0.*dth0dx(1)/(Pe*D);
% c3 = c1 + ddth0dx(1)/Pe;
% c4 = c2-th1hat0(1);
% c3 = th1hat0(1)+c2;
% c3=-1; c4=-1;
S = 0;

%S = c1 + c3.*Pintildeatt + exp(-X).*th1hat0.*(sin(t)+(1/Pe)*cos(t))+ddth0dx(1)/Pe;
yt = Fw + S;
%yt = yt'; 

end