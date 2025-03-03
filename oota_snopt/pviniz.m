function pviniz(tof, r1, v1, dv1, dv2)

% primer vector initialization

% input

%  tof = time-of-flight between impulses
%  r1  = initial position vector (kilometers)
%  v1  = initial velocity vector (kilometers/second)
%  dv1 = initial delta-v vector (kilometers/second)
%  dv2 = final delta-v vector (kilometers/second)

% output via global

%  pvi  = initial primer vector
%  pvdi = initial primer vector derivative

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu pvi pvdi

% compute primer vector at first impulse
  
dv1m = norm(dv1);

pvi = dv1 / dv1m;

% compute primer vector at second impulse
  
dv2m = norm(dv2);

pvf = dv2 / dv2m;
 
% compute state transition matrix

[~, ~, stm] = stm2(mu, tof, r1, v1);

% extract submatrices of state transition matrix

stm11(1:3, 1:3) = stm(1:3, 1:3);

stm12(1:3, 1:3) = stm(1:3, 4:6);

% compute initial value of primer derivative vector

pvdi = stm12 \ (pvf' - stm11 * pvi');

