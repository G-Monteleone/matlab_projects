function [teci2op, top2eci] = tmatrices(inc, raan) 

% orbit plane to/from eci transformation matrices

% input

%  inc  = orbital inclinations (radians)
%  raan = raans (radians) 

% output

%  teci2op = eci-to-orbit plane transformation matrix
%  top2eci = orbit-plane-to-eci transformation matrix

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eci u vector

ueci(1) = cos(raan(2));
ueci(2) = sin(raan(2));
ueci(3) = 0.0;

% eci w vector - first impulse

weci1(1) = sin(raan(1)) * sin(inc(1));
weci1(2) = -cos(raan(1)) * sin(inc(1));
weci1(3) = cos(inc(1));

% eci w vector - second impulse

weci2(1) = sin(raan(2)) * sin(inc(2));
weci2(2) = -cos(raan(2)) * sin(inc(2));
weci2(3) = cos(inc(2));

% eci n vector

xneci = cross(weci2, weci1);

xnecim = sqrt(dot(xneci, xneci));

% eci n unit vector
      
if (xnecim == 0.0)

   xneci(1) = 1.0;
   xneci(2) = 0.0;
   xneci(3) = 0.0;

else

   xneci = xneci / xnecim;

end

xndotu = dot(xneci, ueci);

% argument of hinge line

xphi = acos(xndotu);

if (xneci(3) < 0.0)

   xphi = -xphi;

end

% reference coordinate system-to-eci transformation matrix

top2eci(1, 1) = cos(raan(2)) * cos(xphi) ...
   - sin(raan(2)) * cos(inc(2)) * sin(xphi);

top2eci(1, 2) = -cos(raan(2)) * sin(xphi) ...
   - sin(raan(2)) * cos(inc(2)) * cos(xphi);

top2eci(1, 3) = sin(raan(2)) * sin(inc(2));

top2eci(2, 1) = sin(raan(2)) * cos(xphi) ...
   + cos(raan(2)) * cos(inc(2)) * sin(xphi);

top2eci(2, 2) = -sin(raan(2)) * sin(xphi) ...
   + cos(raan(2)) * cos(inc(2)) * cos(xphi);

top2eci(2, 3) = -cos(raan(2)) * sin(inc(2));

top2eci(3, 1) = sin(inc(2)) * sin(xphi);
top2eci(3, 2) = sin(inc(2)) * cos(xphi);
top2eci(3, 3) = cos(inc(2));

% eci-to-reference coordinate system transformation matrix

teci2op = top2eci';