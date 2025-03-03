% ipto_matlab_snopt.m        January 11, 2022

% two impulse ballistic interplanetary trajectory optimization

% JPL DE421 ephemeris, SNOPT and MICE routines

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global iephem km ephname eq2000

global smu ip1 ip2 jdtdb0 jdtdb1 jdtdb2

global dtr otype imcon revmax aunit boev

global dv_depart dv_arrive ri vi rf vf vito

% initialize jpl ephemeris

ephname = 'de421.bsp';

iephem = 1;

km = 1;

% J2000 equatorial-to-ecliptic transformation matrix

eq2000(1, 1) =  1.000000000000000;
eq2000(1, 2) =  4.403600000000000e-007;
eq2000(1, 3) = -1.909190000000000e-007;

eq2000(2, 1) = -4.799660000000000e-007;
eq2000(2, 2) =  0.917482137087000;
eq2000(2, 3) = -0.397776982902000;

eq2000(3, 1) = 0.000000000000000;
eq2000(3, 2) = 0.397776982902000;
eq2000(3, 3) = 0.917482137087000;

% angular conversion factors

rtd = 180.0 / pi;

dtr = pi / 180.0;

atr = dtr / 3600.0;

% gravitational constant of the sun (km^3/sec^2)

smu = 132712440040.944;

% astronomical unit (kilometers)

aunit = 149597870.70;

% define "reference" julian date (1/1/2000)

jdtdb0 = 2451544.5;

% define celestial body name vector

pname = ['Mercury       '; 'Venus         '; 'Earth         '; 'Mars          '; ...
    'Jupiter       '; 'Saturn        '; 'Uranus        '; 'Neptune       '; ...
    'Pluto         '; 'asteroid/comet'];

% begin simulation

clc; home;

fprintf('\nprogram ipto_matlab_snopt\n');

fprintf('\n< interplanetary trajectory optimization >\n\n');

% request departure calendar date

fprintf('\ndeparture conditions - start date\n');

[month, day, year] = getdate;

jdtdb1 = julian(month, day, year);

while(1)
    
    fprintf('\nplease input the departure date search boundary in days\n');
    
    ddays1 = input('? ');
    
    if (ddays1 >= 0.0)
        
        break;
        
    end
    
end

% request arrival calendar date

fprintf('\n\narrival conditions - start date\n');

[month, day, year] = getdate;

jdtdb2 = julian(month, day, year);

% request search information

while(1)
    
    fprintf('\nplease input the arrival date search boundary in days\n');
    
    ddays2 = input('? ');
    
    if (ddays2 >= 0.0)
        
        break;
        
    end
    
end

% request departure and arrival pecl2000s

for i = 1:1:2
    
    fprintf('\n celestial body menu\n');
    
    fprintf('\n  <1>  Mercury');
    fprintf('\n  <2>  Venus');
    fprintf('\n  <3>  Earth');
    fprintf('\n  <4>  Mars');
    fprintf('\n  <5>  Jupiter');
    fprintf('\n  <6>  Saturn');
    fprintf('\n  <7>  Uranus');
    fprintf('\n  <8>  Neptune');
    fprintf('\n  <9>  Pluto');
    fprintf('\n  <10> asteroid/comet');
    
    if (i == 1)
        
        while(1)
            
            fprintf('\n\nplease select the departure celestial body\n');
            
            ip1 = input('? ');
            
            if (ip1 >= 1 && ip1 <= 10)
                
                break;
                
            end
            
        end
        
    end
    
    if (i == 2)
        
        while(1)
            
            fprintf('\n\nplease select the arrival celestial body\n');
            
            ip2 = input('? ');
            
            if (ip2 >= 1 && ip2 <= 10)
                
                break;
                
            end
            
        end
        
    end
    
end

if (ip2 == 10)
    
    [datafile, pathname] = uigetfile('*.dat', 'please input the name of the asteroid/comet data file');
    
    fid = fopen(datafile, 'r');
    
    % read 25 lines of information
    
    for i = 1:1:25
        
        cline1 = fgetl(fid);
        
        switch i
            
            case 10
                
                tl = size(cline1);
                
                ci = strfind(cline1, ',');
                
                % extract month, day and year
                
                month = str2double(cline1(1:ci(1)-1));
                
                day = str2double(cline1(ci(1)+1:ci(2)-1));
                
                year = str2double(cline1(ci(2)+1:tl(2)));
                
            case 13
                
                boev(1) = str2double(cline1);
                
            case 16
                
                boev(2) = str2double(cline1);
                
            case 19
                
                boev(3) = str2double(cline1);
                
            case 22
                
                boev(4) = str2double(cline1);
                
            case 25
                
                boev(5) = str2double(cline1);
                
        end
        
    end
    
    % julian date of perihelion passage
    
    boev(6) = julian(month, day, year);
    
    fclose(fid);
end

% less than one rev transfer

revmax = 0;

% mission constraints option

while(1)
    
    fprintf('\nwould you like to enforce mission constraints (y = yes, n = no)\n');
    
    slct = input('? ', 's');
    
    if (slct == 'y' || slct == 'n')
        
        break;
        
    end
    
end

imcon = 0;

if (slct == 'y')
    
    imcon = 1;
    
    while(1)
        
        fprintf('\nplease input the lower bound for departure C3 (kilometers^2/second^2)\n');
        
        flow_c3 = input('? ');
        
        if (flow_c3 > 0.0)
            
            break;
            
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the upper bound for departure C3 (kilometers^2/second^2)\n');
        
        fupp_c3 = input('? ');
        
        if (fupp_c3 > flow_c3)
            
            break;
            
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the lower bound for departure DLA (degrees)\n');
        
        flow_dla = input('? ');
        
        if (flow_dla >= -90.0)
            
            break;
            
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the upper bound for departure DLA (degrees)\n');
        
        fupp_dla = input('? ');
        
        if (fupp_dla <= +90.0)
            
            break;
            
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the lower bound for time-of-flight (days)\n');
        
        flow_tof = input('? ');
        
        if (flow_tof > 0.0)
            
            break;
            
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the upper bound for time-of-flight (days)\n');
        
        fupp_tof = input('? ');
        
        if (fupp_tof > flow_tof)
            
            break;
            
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the lower bound for arrival v-infinity (kilometers/second)\n');
        
        flow_vinf = input('? ');
        
        if (flow_vinf > 0.0)
            
            break;
            
        end
        
    end
    
    while(1)
        
        fprintf('\nplease input the upper bound for arrival v-infinity (kilometers/second)\n');
        
        fupp_vinf = input('? ');
        
        if (fupp_vinf > flow_vinf)
            
            break;
            
        end
        
    end
    
end

% request type of optimization

while(1)
    
    fprintf('\n      optimization menu\n');
    
    fprintf('\n <1> minimize departure delta-v\n');
    
    fprintf('\n <2> minimize arrival delta-v\n');
    
    fprintf('\n <3> minimize total delta-v\n');
    
    fprintf('\n <4> no optimization\n');
    
    fprintf('\n selection (1, 2, 3 or 4)\n');
    
    otype = input('? ');
    
    if (otype == 1 || otype == 2 || otype == 3 || otype == 4)
        
        break;
        
    end
    
end

if (otype < 4)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % find optimal solution %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    xg = zeros(2, 1);
    
    xg(1) = jdtdb1 - jdtdb0;
    
    xg(2) = jdtdb2 - jdtdb0;
    
    % bounds on control variables
    
    xlwr = zeros(2, 1);    
    xupr = zeros(2, 1);
    
    xlwr(1) = xg(1) - ddays1;
    xupr(1) = xg(1) + ddays1;
    
    xlwr(2) = xg(2) - ddays2;
    xupr(2) = xg(2) + ddays2;
    
    if (imcon == 1)
        
        flow = zeros(5, 1);
        fupp = zeros(5, 1);
        
        % bounds on nonlinear constraints
        
        flow(2) = flow_c3;
        fupp(2) = fupp_c3;
        
        flow(3) = dtr * flow_dla;
        fupp(3) = dtr * fupp_dla;
        
        flow(4) = flow_tof;
        fupp(4) = fupp_tof;
        
        flow(5) = flow_vinf;
        fupp(5) = fupp_vinf;
        
        fmul = zeros(5, 1);

        fstate = zeros(5, 1);
    
    else
        
        fmul = zeros(1, 1);

        fstate = zeros(1, 1);   
        
    end
   
    % bounds on objective function
    
    flow(1) = 0.0;
    fupp(1) = +Inf;
    
    xmul = zeros(2, 1);

    xstate = zeros(2, 1);
    
    % find optimum
    
    snscreen on;
    
    [x, ~, inform, xmul, fmul] = snopt(xg, xlwr, xupr, xmul, xstate, ...
        flow, fupp, fmul, fstate, 'iptofunc');
        
    if (imcon == 1 && inform ~= 1)
        
        fprintf('\n\ncheck solution!!\n');
        
        fprintf('\n\nall mission constraints may not be satisfied\n');
        
    end
    
    % solution tdb julian dates
    
    jdtdb1 = x(1) + jdtdb0;
    
    jdtdb2 = x(2) + jdtdb0;
    
    % transfer time (days)
    
    taud = jdtdb2 - jdtdb1;
    
    % evaluate current solution
    
    [f, g] = iptofunc (x);
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % no optimization - solve TPBVP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    taud = jdtdb2 - jdtdb1;
    
    tof = taud * 86400.0;
    
    % compute initial state vector
    
    [ri, vi] = becl2000(ip1, jdtdb1);
    
    % compute final state vector
    
    [rf, vf] = becl2000(ip2, jdtdb2);
    
    % solve Lambert's problem
    
    sv1(1:3) = ri;
    
    sv1(4:6) = vi;
    
    sv2(1:3) = rf;
    
    sv2(4:6) = vf;
    
    [vito, vfto] = glambert(smu, sv1, sv2, tof, revmax);
    
    % calculate departure delta-v
    
    dv_depart(1) = vito(1) - vi(1);
    dv_depart(2) = vito(2) - vi(2);
    dv_depart(3) = vito(3) - vi(3);
    
    dvm_depart = norm(dv_depart);
    
    % calculate arrival delta-v
    
    dv_arrive(1) = vf(1) - vfto(1);
    dv_arrive(2) = vf(2) - vfto(2);
    dv_arrive(3) = vf(3) - vfto(3);
    
    dvm_arrive = norm(dv_arrive);
    
end

% convert solution julian dates to calendar dates and tdb

[cdstr1, utstr1] = jd2str(jdtdb1);

[cdstr2, utstr2] = jd2str(jdtdb2);

fprintf('\nDEPARTURE CONDITIONS');
fprintf('\n====================\n');

fprintf('\ndeparture celestial body     ');

disp(pname(ip1, 1:14));

fprintf('\ndeparture calendar date      ');

disp(cdstr1);

fprintf('\ndeparture TDB time           ');

disp(utstr1);

fprintf('\ndeparture julian date        %12.6f\n', jdtdb1);

if (ip1 == 3)
    
    fprintf('\ndeparture delta-v and energy requirements');
    fprintf('\n(Earth mean equator and equinox of J2000)');
    fprintf('\n-----------------------------------------\n');
    
else
    
    fprintf('\ndeparture delta-v and energy requirements');
    fprintf('\n(mean ecliptic and equinox of J2000)');
    fprintf('\n------------------------------------\n');
    
end

if (ip1 == 3)
    
    % compute departure velocity vector in equatorial frame
    
    dveq1 = eq2000 * dv_depart;
    
    dveqm1 = norm(dveq1);
    
    % compute orientation of the departure hyperbola
    
    decl1 = 0.5 * pi - acos(dveq1(3) / dveqm1);
    
    rasc1 = atan3(dveq1(2), dveq1(1));
    
else
    
    decl1 = 0.5 * pi - acos(dv_depart(3) / norm(dv_depart));
    
    rasc1 = atan3(dv_depart(2), dv_depart(1));
    
end

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dveq1(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dveq1(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dveq1(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second', 1000.0 * dveqm1);

fprintf('\n\nenergy                      %12.6f  kilometers^2/second^2\n', dveqm1 * dveqm1);

fprintf('\nasymptote right ascension   %12.6f  degrees\n', rtd * rasc1);

fprintf('\nasymptote declination       %12.6f  degrees\n', rtd * decl1);

fprintf('\nheliocentric orbital elements and state vector of departure celestial body');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

[rplanet, vplanet] = p2000_ecl(ip1, jdtdb1);

oev = eci2orb1(smu, rplanet, vplanet);

oeprint1(smu, oev, 3);

svprint(rplanet, vplanet);

% print orbital elements and state vector of the initial orbit

fprintf('spacecraft heliocentric orbital elements and state vector prior to the first maneuver');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev = eci2orb1(smu, ri, vi);

oeprint1(smu, oev, 3);

svprint(ri, vi);

% print orbital elements and state vector of the transfer orbit

fprintf('spacecraft heliocentric orbital elements and state vector after the first maneuver');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev = eci2orb1(smu, ri, vito');

oeprint1(smu, oev, 3);

svprint(ri, vito);

rito = ri;

fprintf('spacecraft heliocentric orbital elements and state vector prior to the second maneuver');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

tau = 86400.0 * (jdtdb2 - jdtdb1);

[rft, vft] = twobody2(smu, tau, ri, vito');

oev = eci2orb1(smu, rft, vft);

oeprint1(smu, oev, 3);

svprint(rft, vft);

% print orbital elements and state vector of the final orbit

fprintf('spacecraft heliocentric orbital elements and state vector after the second maneuver');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev = eci2orb1 (smu, rft, vft + dv_arrive);

oeprint1(smu, oev, 3);

svprint(rft, vft + dv_arrive');

fprintf('transfer time               %12.6f  days \n ', taud);

% compute and display transfer angle (degrees)

tangle = transfer_angle(rito, rft);

fprintf('\ntransfer angle              %12.6f  degrees \n ', rtd * tangle);

fprintf('\nARRIVAL CONDITIONS');
fprintf('\n==================\n');

% orbital elements and state vector of the arrival body

[rfb, vfb] = p2000_ecl(ip2, jdtdb2);

fprintf('\narrival celestial body       ');

disp(pname(ip2, 1:14));

fprintf('\narrival calendar date        ');

disp(cdstr2);

fprintf('\narrival TDB time             ');

disp(utstr2);

fprintf('\narrival julian date          %12.6f\n', jdtdb2');

fprintf('\nheliocentric orbital elements and state vector of arrival body');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n------------------------------------\n');

oev = eci2orb1(smu, rfb, vfb);

oeprint1(smu, oev, 3);

svprint(rfb, vfb);

% compute arrival v-infinity velocity vector in equatorial frame

dveq2 = eq2000 * (-dv_arrive);

dveqm2 = norm(dveq2);

dveq2prt = dveq2;

if (ip2 == 4)
    
    % Mars mean equator and IAU node of epoch
    
    % compute arrival v-infinity velocity vector in equatorial frame
    
    dveq2 = eq2000 * (-dv_arrive);
    
    dveqm2 = norm(dveq2);
    
    dveq2prt = dveq2;
    
    tmatrix = mme2000(jdtdb2);
    
    dvwrk = dveq2;
    
    dveq2 = tmatrix * dvwrk;
    
end

% compute orientation of the arrival hyperbola

decl2 = 0.5 * pi - acos(dveq2(3) / dveqm2);

rasc2 = atan3(dveq2(2), dveq2(1));

fprintf('arrival delta-v and energy requirements');
fprintf('\n(mean ecliptic and equinox of J2000)');
fprintf('\n-----------------------------------\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dveq2prt(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dveq2prt(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dveq2prt(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second', 1000.0 * dveqm2);

fprintf('\n\nenergy                      %12.6f  kilometers^2/second^2\n', dveqm2 * dveqm2);

if (ip2 == 4)
    
    fprintf('\nMars-mean-equator and IAU node of epoch');
    fprintf('\n---------------------------------------\n');
    
end

fprintf('\nasymptote right ascension   %12.6f  degrees\n', rtd * rasc2);

fprintf('\nasymptote declination       %12.6f  degrees\n', rtd * decl2);

fprintf('\n\ntotal delta-v               %12.6f  meters/second', 1000.0 * (dveqm1 + dveqm2));

fprintf('\n\ntotal energy                %12.6f  kilometers^2/second^2\n', ...
    dveqm1 * dveqm1 + dveqm2 * dveqm2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trajectory graphics options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while(1)
    
    fprintf('\n\nwould you like to plot this trajectory (y = yes, n = no)\n');
    
    slct = input('? ', 's');
    
    if (slct == 'y' || slct == 'n')
        
        break;
        
    end
    
end

if (slct == 'y')
    
    while(1)
        
        fprintf('\nplease input the plot step size (days)\n');
        
        deltat = input('? ');
        
        if (deltat > 0.0)
            
            break;
            
        end
        
    end
    
    % compute orbital periods
    
    [r1, v1] = p2000_ecl(ip1, jdtdb1);
    
    oev1 = eci2orb1(smu, r1, v1);
    
    period1 = 2.0 * pi * oev1(1) * sqrt(oev1(1) / smu) / 86400.0;
    
    [r2, v2] = p2000_ecl(ip2, jdtdb1);
    
    oev2 = eci2orb1(smu, r2, v2);
    
    period2 = 2.0 * pi * oev2(1) * sqrt(oev2(1) / smu) / 86400.0;
    
    xve = oev1(1) / aunit;
    
    if (oev2(1) > oev1(1))
        
        xve = oev2(1) / aunit;
        
    end
    
    % determine number of data points to plot
    
    npts1 = fix(period1 / deltat);
    
    x1 = zeros(npts1 + 1, 1);
    
    y1 = zeros(npts1 + 1, 1);
    
    npts2 = fix(period2 / deltat);
    
    x2 = zeros(npts2 + 1, 1);
    
    y2 = zeros(npts2 + 1, 1);
    
    npts3 = fix((jdtdb2 - jdtdb1) / deltat);
    
    x3 = zeros(npts3 + 1, 1);
    
    y3 = zeros(npts3 + 1, 1);
    
    [rti, vti] = orb2eci(smu, oev);
    
    % create departure orbit data points
    
    for i = 0:1:npts1
        
        jdate = jdtdb1 + i * deltat;
        
        [r1, ~] = p2000_ecl(ip1, jdate);
        
        x1(i + 1) = r1(1) / aunit;
        
        y1(i + 1) = r1(2) / aunit;
        
    end
    
    % compute last data point
    
    [r1, v1] = p2000_ecl(ip1, jdtdb1 + period1);
    
    x1(npts1 + 1) = r1(1) / aunit;
    
    y1(npts1 + 1) = r1(2) / aunit;
    
    % create arrival orbit data points
    
    for i = 0:1:npts2
        
        jdate = jdtdb1 + i * deltat;
        
        [r2, ~] = p2000_ecl(ip2, jdate);
        
        x2(i + 1) = r2(1) / aunit;
        
        y2(i + 1) = r2(2) / aunit;
        
    end
    
    % compute last data point
    
    [r2, ~] = p2000_ecl(ip2, jdtdb1 + period2);
    
    x2(npts2 + 1) = r2(1) / aunit;
    
    y2(npts2 + 1) = r2(2) / aunit;
    
    % create transfer orbit data points
    
    for i = 0:1:npts3
        
        tau = 86400.0 * i * deltat;
        
        [rft, ~] = twobody2(smu, tau, ri, vito');
        
        x3(i + 1) = rft(1) / aunit;
        
        y3(i + 1) = rft(2) / aunit;
        
    end
    
    % compute last data point
    
    tau = 86400.0 * (jdtdb2 - jdtdb1);
    
    [rft, vft] = twobody2(smu, tau, ri, vito');
    
    x3(npts3 + 1) = rft(1) / aunit;
    
    y3(npts3 + 1) = rft(2) / aunit;
    
    % plot orbits and transfer trajectory
    
    figure(1);
    
    hold on;
    
    plot(x1, y1, '.b');
    
    plot(x1, y1, '-b');
    
    plot(x2, y2, '.g');
    
    plot(x2, y2, '-g');
    
    plot(x3, y3, '.r');
    
    plot(x3, y3, '-r');
    
    % plot and label vernal equinox direction
    
    line ([0.05, 1.1 * xve], [0, 0], 'Color', 'black');
    
    text(1.15 * xve, 0, '\Upsilon');
    
    % label launch and arrival locations
    
    [r2, ~] = p2000_ecl(ip2, jdtdb1);
    
    plot(r2(1) / aunit, r2(2) / aunit, '*r');
    
    text(r2(1) / aunit + 0.05, r2(2) / aunit + 0.05, 'L');
    
    [r2, v2] = p2000_ecl(ip1, jdtdb2);
    
    plot(r2(1) / aunit, r2(2) / aunit, '*r');
    
    text(r2(1) / aunit + 0.05, r2(2) / aunit + 0.05, 'A');
    
    plot(x3(1), y3(1), '*r');
    
    text(x3(1) + 0.05, y3(1) + 0.05, 'L');
    
    plot(x3(npts3 + 1), y3(npts3 + 1), '*r');
    
    text(x3(npts3 + 1) + 0.05, y3(npts3 + 1) + 0.05, 'A');
    
    % label launch and arrival dates
    
    text(0.85 * xve, -xve + 0.8, 'Launch ', 'FontSize', 8);
    
    text(0.875 * xve, -xve + 0.7, pname(ip1, 1:14), 'FontSize', 8);
    
    text(0.875 * xve, -xve + 0.6, cdstr1, 'FontSize', 8);
    
    text(0.85 * xve, -xve + 0.4, 'Arrival', 'FontSize', 8);
    
    text(0.875 * xve, -xve + 0.3, pname(ip2, 1:14), 'FontSize', 8);
    
    text(0.875 * xve, -xve + 0.2, cdstr2, 'FontSize', 8);
    
    % label plot and axes
    
    xlabel('X coordinate (AU)', 'FontSize', 12);
    
    ylabel('Y coordinate (AU)', 'FontSize', 12);
    
    title('Interplanetary Trajectory Optimization', 'FontSize', 16);
    
    % plot aphelion and perihelion of departure planet
    
    oev1(6) = 0.0;
    
    [r, ~] = orb2eci(smu, oev1);
    
    plot(r(1) / aunit, r(2) / aunit, 'sb');
    
    oev1(6) = pi;
    
    [r, ~] = orb2eci(smu, oev1);
    
    plot(r(1) / aunit, r(2) / aunit, 'ob');
    
    % plot aphelion and perihelion of arrival planet
    
    oev2(6) = 0.0;
    
    [r, ~] = orb2eci(smu, oev2);
    
    plot(r(1) / aunit, r(2) / aunit, 'sg');
    
    oev2(6) = pi;
    
    [r, v] = orb2eci(smu, oev2);
    
    plot(r(1) / aunit, r(2) / aunit, 'og');
    
    if (ip1 == 3)
        
        % plot line of nodes (Earth is departure planet)
        
        oev2(6) = -oev2(4);
        
        [r, ~] = orb2eci(smu, oev2);
        
        x4(1) = r(1) / aunit;
        
        y4(1) = r(2) / aunit;
        
        oev2(6) = oev2(6) + pi;
        
        [r, v] = orb2eci(smu, oev2);
        
        x4(2) = r(1) / aunit;
        
        y4(2) = r(2) / aunit;
        
        plot(x4, y4, ':g');
        
    end
    
    % plot sun
    
    plot(0, 0, 'hy', 'MarkerSize', 10);
    
    axis equal;
    
    zoom on;
    
    % the next line creates a color eps graphics file with tiff preview
    
    print -depsc -tiff -r300 ipto_matlab.eps
    
end

while(1)
    
    fprintf('\nwould you like to plot the primer characteristics (y = yes, n = no)\n');
    
    slct = input('? ', 's');
    
    if (slct == 'y' || slct == 'n')
        
        break;
        
    end
    
end

if (slct == 'y')
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % primer vector analysis
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % perform primer vector initialization
    
    tof = 86400.0 * (jdtdb2 - jdtdb1);
    
    [pvi, pvdi] = pv_iniz(smu, tof, rito, vito', dv_depart, dv_arrive);
    
    % number of graphic data points
    
    npts = 300;
    
    % plot behavior of primer vector magnitude
    
    dt = tof / npts;
    
    xp1 = zeros(npts + 1, 1);
    
    yp1 = zeros(npts + 1, 1);
    
    yp2 = zeros(npts + 1, 1);
    
    for i = 1:1:npts + 1
        
        t = (i - 1) * dt;
        
        if (t == 0.0)
            
            % initial value of primer magnitude and derivative
            
            pvm = norm(pvi);
            
            pvdm = dot(pvi, pvdi) / pvm;
            
        else
            
            % primer vector and derivative magnitudes at time t
            
            [pvm, pvdm] = pvector(smu, rito, vito', pvi, pvdi, t);
            
        end
        
        % load data array
        
        xp1(i) = t / 86400.0;
        
        yp1(i) = pvm;
        
        yp2(i) = pvdm;
        
    end
    
    figure(2);
    
    hold on;
    
    plot(xp1, yp1, '-r', 'LineWidth', 1.5);
    
    plot(xp1, yp1, '.r');
    
    title('Primer Vector Analysis', 'FontSize', 16);
    
    xlabel('simulation time (days)', 'FontSize', 12);
    
    ylabel('primer vector magnitude', 'FontSize', 12);
    
    grid;
    
    print -depsc -tiff -r300 primer.eps;
    
    % plot behavior of magnitude of primer derivative
    
    figure(3);
    
    hold on;
    
    plot(xp1, yp2, '-r');
    
    plot(xp1, yp2, '.r');
    
    title('Primer Vector Analysis', 'FontSize', 16);
    
    xlabel('simulation time (days)', 'FontSize', 12);
    
    ylabel('primer derivative magnitude', 'FontSize', 12);
    
    grid;
    
    print -depsc -tiff -r300 primer_der.eps;
    
end

fprintf('\n\n');

% CSPICE_KCLEAR clears the KEEPER system: unload all kernels, clears
% the kernel pool, and re-initialize the system.

cspice_kclear;


