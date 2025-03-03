function rrd = jpleph_mice(et, ntarg, ncent)

% reads the jpl planetary ephemeris and gives the position and velocity
% of the point 'ntarg' with respect to point 'ncent' using MICE routines

% input

%   et    = TDB julian date at which interpolation is wanted

%   ntarg = integer number of 'target' point

%   ncent = integer number of center point

%   the numbering convention for 'ntarg' and 'ncent' is:

%   the numbering convention for 'ntarg' and 'ncent' is:

%        1 = mercury           8 = neptune
%        2 = venus             9 = pluto
%        3 = earth            10 = moon
%        4 = mars             11 = sun
%        5 = jupiter          12 = solar-system barycenter
%        6 = saturn           13 = earth-moon barycenter
%        7 = uranus           14 = nutations (longitude and obliq)
%                             15 = librations, if on ephemeris file

% output

%   rrd = output 6-word array containing position and velocity
%         of point 'ntarg' relative to 'ncent'. the units are
%         determined by the value of km passed via global.

% global

%   iephem  = initialization flag (1 = initialize)
%   ephname = name of ephemeris binary data file (de421.bsp, etc.)
%   km      = state vector units flag (1 = km & km/sec, 0 = au & au/day)
%   aunit   = numerical value of astronomical unit (kilometers)

% NOTE: output in Earth true-of-date system

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global iephem ephname km aunit

if (iephem == 1)
    
    % load binary ephemeris data file
    
    cspice_furnsh(ephname);
    
    % reset initialization flag
    
    iephem = 0;
    
end

% set name of target body

switch ntarg
   
    case (1)
        
        targ = 'mercury';
        
    case (2)
        
        targ = 'venus';
        
    case (3)
        
        targ = 'earth';
        
    case (4)
        
        targ = 'mars';
        
    case(5)
        
        targ = '5';
        
    case (6)
        
        targ = '6';
        
    case (7)
        
        targ = '7';
        
    case (8)
        
        targ = '8';
        
    case (9)
        
        targ = '9';
        
    case (10)
        
        targ = 'moon';
        
    case (11)
        
        targ = 'sun';
        
    case (15)
        
        targ = '15';
        
end

% set name of central body

switch ncent
        
    case (1)
        
        obs = 'mercury';
        
    case (2)
        
        obs = 'venus';
        
    case (3)
        
        obs = 'earth';
        
    case (4)
        
        obs = 'mars';
        
    case (5)
        
        obs = '5';
        
    case (6)
        
        obs = '6';
        
    case (7)
        
        obs = '7';
        
    case (8)
        
        obs = '8';
        
    case (9)
        
        obs = '9';
        
    case (10)
        
        obs = 'moon';
        
    case (11)
        
        obs = 'sun';
        
    case (12)
        
        % solar system barycenter
        
        obs = '0';       
end

% compute time, expressed as TDB seconds past J2000 TDB (2451545.0)

etime = 86400.0 * (et - 2451545.0);

% compute position and velocity vectors in ecliptic j2000 system (no corrections)

starg = mice_spkezr(targ, etime, 'ECLIPJ2000', 'NONE', obs);

% provide output in user-requested units

if (km == 1)
    
    % state is kilometers and kilometers/second
    
    rrd = starg.state;
    
else
    
    % state is aunit and aunit/day
    
    rrd(1:3) = starg.state(1:3) / aunit;
    
    rrd(4:6) = 86400.0 * starg.state(4:6) / aunit;
    
    rrd = rrd';
    
end





