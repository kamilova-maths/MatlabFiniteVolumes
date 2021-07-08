global Pe Bi tha N K Gamma St T L D uf x1 x2 Q Ldim uc c1 c1dim P0 eps

params = 'stand';
switch(params)
    case 'stand'
        rho= 1.8*10^3; %Bergstrom ; 
        g = 10; 
        c= 900; % Fitt and Howell
        Ldim=5; % Temperature Profiles in Soderberg Electrodes
        uc = 10^-5; %Bergstrom approximation
        %uc = 5.4*10^(-3); % Periodic value (where the total sum still corresponds to Bergstrom approx.
        %uc = 0.5/(24*60*60); % actual estimated value
        R0=0.5;
        R1=1; 
        k = 3; 
        Qc=15000; %Taken very vaguely from Temperature profiles in Soderberg electrodes. 
        %Qc = 3.5607e+03;
        %Qc = 2000/L; 
        %mu0=10^10; % given by Bjornar at a reference temperature
        mu0 = 4.5E10; 
        T_in= 20 + 273; 
        T_a = 80 + 273; 
        h = 7; 

        %Defining non-dimensional parameters
        % Peclet number
        Pe = (rho*c*uc*Ldim)/(k);
        epsilon=R1/Ldim;
        St=(rho*g*Ldim^2)/(uc*mu0);
        %St =37;
        P0 = (10000*Ldim)/((R1^2)*uc*mu0);
        %P0 = 0.1;
        %P0 = 1; 
        Bi= ((Ldim^2)*h)/(k*R1); 
        %Bi =100;
        %disp('Remember Bi is 27')
        %Bi = Pe;
        %tha = 0.005; 
        %D = (R0^2)/(R1^2); 
        D = 0.5;
        %D = 0.18;
        %This is the area of the clamps, taken from Temperature profiles ... 
        x1dim = 5-2;
        x2dim = 6.5-2; 
        %x1 = x1dim/Ldim;
        %x2 = x2dim/Ldim;

        x1 = 0.7;
        x2 = 0.9;
        % Q scale is Q0*Qprime
        eps = 1e-4;
        %DeltaT = (Qc*Ldim)/(rho*c*uc);
        DeltaT = 280+273-T_in; 
        Qscale = DeltaT*(rho*c*uc)/(Ldim); 
        Q0 = Qscale*(x2dim-x1dim)/Ldim; % alternatively, Q0 can be whatever we make it

        %Q0 = Qc*(x2dim - x1dim)/Ldim ; 

        % DeltaT =  (Q0*Ldim)/(rho*c*uc); 
        %gammaBar = 0.069; 
        gammaBar = 0.09; 
        Gamma = gammaBar*DeltaT; 
        %Gamma = 23 ;
        %Q = Ldim/(x2dim - x1dim);
        Q = 1/(x2-x1); 
        tha = (T_a- T_in)/DeltaT; 
        %tha=0;
        uf = 1; 
        %uf = 0.5/(uc*24*60*60); 
        %uf = 0.3456/(uc*24*60*60); 
        %uf = 0.6912/(uc*24*60*60); 
        %y0(1+K:2*K) = th0;
        %c1dim = 0.5; 
        c1dim = 1; % Aasgeir, 17/08/2020 - this is the best one to use
        %c1dim = 0.75;
        c1 = c1dim/Ldim; 
        % Pe = 37.8; St = 8.8; P0 =0.7; Bi = 114.3; tha=0.005; D = 0.25; 
        % gamma = 30;  x1 = 5/7; x2 = 6.5/7; Q = 1; uf = 1; 

        % Calculating the initial conditions as a solution of the steady state
        % problem 
        % Discretisation in t
        N = 20000; 
        % Discretisation in x
        K=1500;

        % end of the domain
        %T = 2*pi/0.5;
        L= 3; 
        T = 5;
        
        CFL = 4*(T/(N-1)) / (L/(K-1));
        
        if CFL>=1
            disp('This probably wont work')
        end
        
    case 'diff'
        rho= 1.8*10^3; %Bergstrom ; 
        g = 10; 
        c= 900; % Fitt and Howell
        Ldim=5; % Temperature Profiles in Soderberg Electrodes
        uc = 10^-5; %Bergstrom approximation
        %uc = 5.4*10^(-3); % Periodic value (where the total sum still corresponds to Bergstrom approx.
        %uc = 0.5/(24*60*60); % actual estimated value
        R0=0.5;
        R1=1; 
        k = 3; 
        Qc=15000; %Taken very vaguely from Temperature profiles in Soderberg electrodes. 
        %Qc = 3.5607e+03;
        %Qc = 2000/L; 
        %mu0=10^10; % given by Bjornar at a reference temperature
        mu0 = 4.5E10; 
        T_in= 20 + 273; 
        T_a = 80 + 273; 
        h = 7; 

        %Defining non-dimensional parameters
        % Peclet number
        %Pe = (rho*c*uc*Ldim)/(k);
        Pe = 1;
        epsilon=R1/Ldim;
        %St=(rho*g*Ldim^2)/(uc*mu0);
        St = 1;
        %St =37;
        %P0 = (10000*Ldim)/((R1^2)*uc*mu0);
        %P0 = 0.1;
        P0 = 1; 
        Bi = Pe;  
        %Bi =100;
        %disp('Remember Bi is 27')
        %Bi = Pe;
        %tha = 0.005; 
        D = (R0^2)/(R1^2); 
        

        %This is the area of the clamps, taken from Temperature profiles ... 
        x1dim = 5-2;
        x2dim = 6.5-2; 
        %x1 = x1dim/Ldim;
        %x2 = x2dim/Ldim;
        
        x1 = 0.7;
        x2 = 0.9;
        % Q scale is Q0*Qprime
        eps = 1e-4;
        %DeltaT = (Qc*Ldim)/(rho*c*uc);
        DeltaT = 280+273-T_in; 
        Qscale = DeltaT*(rho*c*uc)/(Ldim); 
        Q0 = Qscale*(x2dim-x1dim)/Ldim; % alternatively, Q0 can be whatever we make it

        %Q0 = Qc*(x2dim - x1dim)/Ldim ; 

        % DeltaT =  (Q0*Ldim)/(rho*c*uc); 
        %gammaBar = 0.069; 
        gammaBar = 0.09; 
        %Gamma = gammaBar*DeltaT; 
        %Gamma = 1; 
        Gamma = 23 ;
        %Q = Ldim/(x2dim - x1dim);
        
        tha = (T_a- T_in)/DeltaT; 
        %tha=0;
        uf = 1; 
        %uf = 0.5/(uc*24*60*60); 
        %uf = 0.3456/(uc*24*60*60); 
        %uf = 0.6912/(uc*24*60*60); 

        %c1dim = 0.5; 
        c1dim = 1; % Aasgeir, 17/08/2020 - this is the best one to use
        %c1dim = 0.75;
        c1 = c1dim/Ldim; 
        % Pe = 37.8; St = 8.8; P0 =0.7; Bi = 114.3; tha=0.005; D = 0.25; 
        % gamma = 30;  x1 = 5/7; x2 = 6.5/7; Q = 1; uf = 1; 

        % Calculating the initial conditions as a solution of the steady state
        % problem 
        % Discretisation in t
        N = 8000; 
        % Discretisation in x
        K=300;

        % end of the domain
        %T = 2*pi/0.5;
        L= 1; 
        b = (x2-x1)./2 + x1; 
        c = x2-x1; 
        Q = (sqrt(2/pi))./(c.*(erf(b/(sqrt(2).*c)) - erf((b-L)./(sqrt(2).*c))));
        Qfun = @(x) Q.*exp(-((x-((x2-x1)./2+x1)).^2)./(2*((x2-x1)).^2));
        %Q = 1/(x2-x1); 
        extra = 0.01; 
        %T = 2*pi+extra;
        %T = 2*0.1728; 
        T = 5;

        % We define non-dimensional day d
        d = 86400*uc/Ldim;
        disp('Bi is equal to Pe, Gamma is what it should, Qfun continuous, everything else is 1, K is 300')
        
end