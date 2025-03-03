% de405_integration_epoch    April 20, 2009

% return the state vector of the Moon as seen from the Earth in
% the J2000 frame at the DE405 integration epoch using MICE routines

% Reference

% E M Standish, JPL IOM 312.F-98-048, August 26, 1998, page 11.

% Orbital Mechanics with Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% Load a leapseconds file and DE405 binary ephemeris

%cspice_furnsh('naif0010.tls');

cspice_furnsh('de405.bsp');

% astronomical unit (kilometers)

aunit = 149597870.691;

format long;

target   = 'Moon';
frame    = 'J2000';
abcorr   = 'NONE';
observer = 'Earth';
         
% number of ephemeris seconds since reference

et = 86400 * (2440400.5 - 2451545.0);

deltat = cspice_deltet(et, 'ET');

% look-up the state for the defined parameters

starg = mice_spkezr(target, et, frame, abcorr, observer);

% display results

disp(' ')

txt = sprintf('The position of the  : %s', target);

disp(txt)

txt = sprintf('As observed from the : %s', observer);

disp(txt)

txt = sprintf('In reference frame   : %s', frame);

disp(txt)

disp(' ')

utc_epoch = cspice_et2utc(et, 'C', 3);

jdutc = cspice_et2utc(et, 'J', 8);

txt = sprintf('At epoch (UTC)     : %s', utc_epoch);

disp(txt)

disp(' ')

txt = sprintf(['R (AU)     : ' ...
    '%18.16f %18.16f %18.16f'], starg.state(1:3)/aunit);

disp(txt)

txt = sprintf(['V (AU/day) : ' ...
    '%18.16f %18.16f %18.16f'], 86400 * starg.state(4:6)/aunit);

disp(txt)

disp(' ')

% unload ephemeris

cspice_unload('de405.bsp')
