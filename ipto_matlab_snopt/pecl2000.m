function [r, v] = pecl2000(ntarg, jdate)

% heliocentric planetary state vector

% earth mean ecliptic and equinox of j2000
% coordinate system (jpl ephemeris)

% input

%  jdate = TDB julian date
%  ntarg = "target" body

% output

%  r = position vector (kilometers)
%  v = velocity vector (kilometers/second)

% NOTE: requires equatorial to ecliptic
%       transformation matrix eq2000 via global

% eq2000 = [+1.000000000000 -0.000000479966  0.000000000000]
%          [+0.000000440360 +0.917482137087 +0.397776982902]
%          [-0.000000190919 -0.397776982902 +0.917482137087]
      
% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global eq2000

% call jpl ephemeris

ncent = 11;

result = jplephem (jdate, ntarg, ncent);

% load position and velocity vectors

rtmp = result(1: 3);

vtmp = result(4: 6);

% convert to ecliptic

r = eq2000 * rtmp;

v = eq2000 * vtmp;



