function f = oota_func2 (x)

% transfer orbit semiparameter objective function

% input

%  x = current x and z values

% output

%  f = objective function = partial(dv)/partial(p)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu slr

global dr r1xr2m dtheta dv1m dv2m 

global tmp1 tmp2 u1 u2

% compute velocity vectors (kilometers/second)

v = sqrt(mu * x) * dr / r1xr2m;
v1 = sqrt(mu / slr(1)) * tmp1;
v2 = sqrt(mu / slr(2)) * tmp2;

z = sqrt(mu / x) * tan(0.5 * dtheta);

% compute delta-v vectors (kilometers/second)

dv1 = (v + z * u1) - v1;
dv2 = v2 - (v - z * u2);

% delta-v magnitudes (kilometers/second)

dv1m = norm(dv1);
      
dv2m = norm(dv2);

t1 = v - z * u1;
t2 = v + z * u2;

% compute dot products

dv1dott1 = dot(dv1, t1);
      
dv2dott2 = dot(dv2, t2);

% objective function = partial(dv)/partial(p)

f = ((dv1dott1/dv1m) - (dv2dott2/dv2m))/(2 * x);
