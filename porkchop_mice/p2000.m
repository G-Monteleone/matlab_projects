function [r, v] = p2000(ncent, ntarg, jdate)

% lunar/planetary state vector

% EME2000 coordinate system (jpl ephemeris)

% input

%  jdate = julian date
%  ncent = "central" body
%  ntarg = "target" body

% output

%  r = position vector
%  v = velocity vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call mice version of jpl ephemeris

result = jpleph_mice (jdate, ntarg, ncent);

% load position and velocity vectors

r = result(1: 3);
   
v = result(4: 6);


   


