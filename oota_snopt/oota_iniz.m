function oota_iniz

% initialization routine

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global rinc

global sma ecc inc argper raan slr

global apeci1 apeci2 teci2op

global aapop1 aapop2 wop1 wop2 ecc1 ecc2

% semiparameters of initial and final orbits

slr(1) = sma(1) * (1.0 - ecc(1) * ecc(1));
      
slr(2) = sma(2) * (1.0 - ecc(2) * ecc(2));

% -----------------------------------------
% compute eci unit angular momentum vectors
% -----------------------------------------

% initial orbit

weci1(1) = sin(raan(1)) * sin(inc(1));
weci1(2) = -cos(raan(1)) * sin(inc(1));
weci1(3) = cos(inc(1));
 
% final orbit

weci2(1) = sin(raan(2)) * sin(inc(2));
weci2(2) = -cos(raan(2)) * sin(inc(2));
weci2(3) = cos(inc(2));

% compute relative inclination (radians)

crinc = dot(weci1, weci2);

rinc = acos(crinc);

% -------------------------------------------------
% compute orbit plane unit angular momentum vectors
% -------------------------------------------------

% initial orbit

wop1(1) = 0;
wop1(2) = -sin(rinc);
wop1(3) = cos(rinc);

% final orbit

wop2(1) = 0;
wop2(2) = 0;
wop2(3) = 1;

% compute eci perigee vectors

apeci1(1) = cos(argper(1)) * cos(raan(1)) ...
   - sin(argper(1)) * sin(raan(1)) * cos(inc(1));

apeci1(2) = cos(argper(1)) * sin(raan(1)) ...
   + sin(argper(1)) * cos(raan(1)) * cos(inc(1));

apeci1(3) = sin(argper(1)) * sin(inc(1));

apeci2(1) = cos(argper(2)) * cos(raan(2)) ...
   - sin(argper(2)) * sin(raan(2)) * cos(inc(2));

apeci2(2) = cos(argper(2)) * sin(raan(2)) ...
   + sin(argper(2)) * cos(raan(2)) * cos(inc(2));

apeci2(3) = sin(argper(2)) * sin(inc(2));

% compute reference coordinate system to/from eci transformation matrices

[teci2op, ~] = tmatrices(inc, raan);

% -------------------------------------------------------
% compute reference coordinate system argument of perigee
% -------------------------------------------------------

apop1 = teci2op * apeci1';

apop2 = teci2op * apeci2';

% first impulse

aapop1 = acos(apop1(1));

if (apop1(3) < 0.0)

   aapop1 = -aapop1;

end

% second impulse

aapop2 = atan3(apop2(2), apop2(1));

% --------------------------------------------------------
% compute reference coordinate system eccentricity vectors
% --------------------------------------------------------

% initial orbit

ecc1(1) = ecc(1) * cos(aapop1);
ecc1(2) = ecc(1) * sin(aapop1) * cos(rinc);
ecc1(3) = ecc(1) * sin(aapop1) * sin(rinc);

% final orbit

ecc2(1) = ecc(2) * cos(aapop2);
ecc2(2) = ecc(2) * sin(aapop2);
ecc2(3) = 0;

