function [pos,vel] = OrbitPropogator(time_JD,planet)

    % GM of sun
    mu = 1.327124e11;                           % [km^3 / s^2]
    
    % unit conversions
    AU2KM = 149597870.691;
    DEG2RAD = pi/180;
    ASEC2DEG = 0.000277778;
    ASEC2RAD = ASEC2DEG * DEG2RAD;
    
    switch planet %Qui potrei pensare di inserire i miei asteroidi con dati acquisiti il 1 gennaio 2000
        case 'Earth'
            
            % J2000 orbital elements
            a = 1.00000011 * AU2KM;             % [km]
            e = 0.01671022;                     % [1]
            i = 0.00005 * DEG2RAD;              % [rad]
            Omega = -11.26064 * DEG2RAD;        % [rad]
            omega = 102.94719 * DEG2RAD;        % [rad]
            L = 100.46435 * DEG2RAD;            % [rad]
            T = 2*pi*sqrt(a^3 / mu);            % [s]

%EC= 1.539046554914033E-02 QR= 9.803298618410823E-01 IN= 1.179971863122956E-02
%  OM= 1.059828342306277E+01 W = 6.373218303624625E+01 Tp=  2451519.103806403000
%  N = 9.927296858629503E-01 MA= 2.570792013466453E+01 TA= 2.648640956461558E+01
%  A = 9.956534316802405E-01 AD= 1.010977001519399E+00 PR= 3.626364811354087E+02

%         case 'Earth' %Sun corretto da horizon con j2000 alle 12
%             
%             % J2000 orbital elements
%             a = 1.00044883 * AU2KM;             % [km]
%             e = 0.01711863;                     % [1]
%             i = 0.000418 * DEG2RAD;             % [rad]
%             Omega = 135.08072 * DEG2RAD;        % [rad]
%             omega = 326.72822 * DEG2RAD;        % [rad]
%             L = 358.61726 * DEG2RAD;            % [rad]
%             T = 2*pi*sqrt(a^3 / mu);            % [s]          

%         case 'Earth' %Sun corretto da horizon con j2000 alle 12
%             
%             % J2000 orbital elements
%             a = 0.99565343 * AU2KM;             % [km]
%             e = 0.01539047;                     % [1]
%             i = 0.01180 * DEG2RAD;              % [rad]
%             Omega = 10.59828 * DEG2RAD;         % [rad]
%             omega = 63.73218 * DEG2RAD;         % [rad]
%             L = 25.70792 * DEG2RAD;             % [rad]
%             T = 2*pi*sqrt(a^3 / mu);            % [s]            

            
            % 1st order fit (via centennial rates)
            da = -0.00000005 * AU2KM;           % [km / JC]
            de = -0.00003804;                   % [1 / JC]
            di = -46.94 * ASEC2RAD;             % [rad / JC]
            dOmega = -18228.25 * ASEC2RAD;      % [rad / JC]
            domega = 1198.28 * ASEC2RAD;        % [rad / JC]
            dL = 129597740.63 * ASEC2RAD;       % [rad / JC]
            
        case 'Mars'
            
            a = 1.52366231 * AU2KM;
            e = 0.09341233;                      % eccentricity
            i = 1.85061 * DEG2RAD;               % inclination
            Omega = 49.57854 * DEG2RAD;          % RAAN
            omega = 336.04084 * DEG2RAD;         % argument of the perigee
            L = 355.45332 * DEG2RAD;             % mean longitude at J2000
            T = 2*pi*sqrt(a^3 / mu);

        case '2006RH120' %3403148 id asteroide su SPICE - DV Da 10 a 15

            a = 0.95317871 * AU2KM;
            e = 0.05071946;                      % eccentricity
            i = 0.56718 * DEG2RAD;               % inclination
            Omega = 311.72627 * DEG2RAD;         % RAAN
            omega = 127.82117 * DEG2RAD;         % argument of the perigee
            L = 207.74711 * DEG2RAD;             % mean longitude at J2000
            T = 2*pi*sqrt(a^3 / mu);

        case '2003YN107' %3170221 id asteroide su SPICE - DV


            %Caso con centro nel sole e non nel SSB
            a = 0.99775910 * AU2KM;
            e = 0.04079613;                      % eccentricity           
            i = 4.27213 * DEG2RAD;               % inclination
            Omega = 276.43156 * DEG2RAD;         % RAAN
            omega = 156.86362 * DEG2RAD;         % argument of the perigee
            L = 27.08720 * DEG2RAD;              % mean longitude at J2000
            T = 2*pi*sqrt(a^3 / mu);

%             a = 0.99356545 * AU2KM;
%             e = 0.04246728;                      % eccentricity
%             i = 4.26927 * DEG2RAD;               % inclination
%             Omega = 277.03770 * DEG2RAD;         % RAAN
%             omega = 145.53831 * DEG2RAD;         % argument of the perigee
%             L = 37.90059 * DEG2RAD;              % mean longitude at J2000
%             T = 2*pi*sqrt(a^3 / mu);            

            

        case '2006JY26' %3332535 id asteroide su SPICE - DV

            a = 0.97959182 * AU2KM;
            e = 0.08706775;                      % eccentricity
            i = 1.98723 * DEG2RAD;               % inclination
            Omega = 47.95470 * DEG2RAD;          % RAAN
            omega = 296.33563 * DEG2RAD;         % argument of the perigee
            L = 57.77792 * DEG2RAD;             % mean longitude at J2000
            T = 2*pi*sqrt(a^3 / mu);

        case '2013BS45' %3625129 id asteroide su SPICE - DV

            a = 1.00845959 * AU2KM;
            e = 0.08474136;                      % eccentricity
            i = 0.93187 * DEG2RAD;               % inclination
            Omega = 98.48003 * DEG2RAD;          % RAAN
            omega = 127.71404 * DEG2RAD;         % argument of the perigee
            L = 299.01174 * DEG2RAD;             % mean longitude at J2000
            T = 2*pi*sqrt(a^3 / mu);

        case '2012FC71' %3602157 id asteroide su SPICE - DV

            a = 0.98758246 * AU2KM;
            e = 0.08814198;                      % eccentricity
            i = 4.97015 * DEG2RAD;               % inclination
            Omega = 39.62314 * DEG2RAD;          % RAAN
            omega = 348.40142 * DEG2RAD;         % argument of the perigee
            L = 352.77063 * DEG2RAD;             % mean longitude at J2000
            T = 2*pi*sqrt(a^3 / mu);

        case '2014EK24' %459872 id asteroide su SPICE - DV

            a = 1.00892236 * AU2KM;
            e = 0.07552598;                      % eccentricity
            i = 4.64463 * DEG2RAD;               % inclination
            Omega = 345.08068 * DEG2RAD;         % RAAN
            omega = 66.00314 * DEG2RAD;          % argument of the perigee
            L = 95.98361 * DEG2RAD;              % mean longitude at J2000
            T = 2*pi*sqrt(a^3 / mu);

    end 
    
    % mean anomaly at J2000
    M_J2000 = L - Omega - omega;
    
    time_J2000 = time_JD - 2451545.0;
    time_sec = time_J2000 * 86400;
     
    Me = M_J2000 + (2*pi/T) * time_sec;
    
    
    E = kepler(e,Me,1e-10);
    theta = 2*atan(tan(E/2) * sqrt((1 + e)/(1 - e)));   % rad
    
    % calculate specific angular momentum
    h = sqrt(a*mu*(1 - e^2));
    
    % position and velocity in perifocal coord. frame
    r_PF = h^2./(mu*(1 + e*cos(theta))) ...
        .* [cos(theta); sin(theta); zeros(size(theta))];
    v_PF = mu/h * [-sin(theta); e + cos(theta); zeros(size(theta))];

    % transform perifocal to heliocentric equatorial system
    % R1 = EulerRotation('z',-omega);
    % R2 = EulerRotation('x',-i);
    % R3 = EulerRotation('z',-Omega);

    R1 = [cos(Omega) sin(Omega) 0;
        -sin(Omega) cos(Omega) 0;
        0 0 1]';

    R2 = [1 0 0;
        0 cos(i) sin(i);
        0 -sin(i) cos(i)]';

    R3 = [cos(omega) sin(omega) 0;
        -sin(omega) cos(omega) 0;
        0 0 1]';

    r_HEC = R1*R2*R3*r_PF;
    v_HEC = R1*R2*R3*v_PF;

    pos = r_HEC';
    vel = v_HEC';
end