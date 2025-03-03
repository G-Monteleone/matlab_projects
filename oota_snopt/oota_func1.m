function f = oota_func1 (x)

% delta-v objective function

% input

%  x(1) = current phi 1 (radians)
%  x(2) = current phi 2 (radians)

% output

%  f = current total delta-v (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global rinc slr ecc dr aapop1 aapop2 ecc1 ecc2

global r1 v1 r2 v2 vt1 vt2 dv1 dv2

global u1 u2 wop1 wop2 tmp1 tmp2

global mu r1xr2m dtheta

% load current values of control variables

phi1 = x(1);
      
phi2 = x(2);

% u vector at first impulse

u1(1) = cos(phi1);
u1(2) = sin(phi1) * cos(rinc);
u1(3) = sin(phi1) * sin(rinc);

% u vector at second impulse

u2(1) = cos(phi2);
u2(2) = sin(phi2);
u2(3) = 0.0;

xtmp1 = slr(1) / (1.0 + ecc(1) * cos(phi1 - aapop1));

xtmp2 = slr(2) / (1.0 + ecc(2) * cos(phi2 - aapop2));

% position vectors at first and second impulse (kilometers)

r1 = xtmp1 * u1;

r2 = xtmp2 * u2;

% geocentric distance at first and second impulses (kilometers)

r1m = norm(r1);
      
r2m = norm(r2);

% transfer orbit true anomaly increment (radians)

dtheta = acos(dot(u1, u2));

dr = r2 - r1;

% perform utility calculations

r1xr2 = cross(r1, r2);

r1xr2m = norm(r1xr2);

tmp1 = cross(wop1, (ecc1 + u1));

tmp2 = cross(wop2, (ecc2 + u2));

r1dotr2 = dot(r1, r2);

% compute transfer orbit semiparameter bounds (kilometers)

pmin = (r1m * r2m - r1dotr2) ...
   / (r1m + r2m + sqrt(2 * (r1m * r2m + r1dotr2)));

pmax = (r1m * r2m - r1dotr2) ...
   / (r1m + r2m - sqrt(2 * (r1m * r2m + r1dotr2)));
    
% compute transfer orbit semiparameter (kilometers)

f1 = oota_func2(pmin);

f2 = oota_func2(pmax);

if (f1 * f2 < 0.0)
    
   rtol = 1.0e-8;

   [xroot, ~] = brent('oota_func2', pmin, pmax, rtol);
   
else
    
   % invalid bracket - set delta-v to large value
   
   f = 1e99;
   
   return;
   
end

% load semiparameter solution

slrt = xroot;

% compute velocity vectors at first and second impulse (kilometers/second)

v = sqrt(mu * slrt) * dr / r1xr2m;

v1 = sqrt(mu / slr(1)) * tmp1;

v2 = sqrt(mu / slr(2)) * tmp2;

% compute delta-v vectors (kilometers/second)

z = sqrt(mu / slrt) * tan(0.5 * dtheta);

dv1 = (v + z * u1) - v1;

dv2 = v2 - (v - z * u2);

% compute delta-v magnitudes (kilometers/second)

dv1m = norm(dv1);
      
dv2m = norm(dv2);

% transfer orbit velocity vectors (kilometers/second)

vt1 = v + z * u1;

vt2 = v - z * u2;

% compute objective function - total delta-v (kilometers/second)

f = dv1m + dv2m;

