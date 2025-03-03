function [pvm, pvdm] = pvector(mu, ri, vi, pvi, pvdi, t)

% primer vector and derivative magnitudes

% input

%  mu   = gravitational constant (kilometers^3/second^2)
%  ri   = current position vector (kilometers)
%  vi   = current velocity vector (kilometers/second)
%  pvi  = initial primer vector
%  pvdi = initial primer vector derivative
%  t    = current time (seconds)

% output

%  pvm  = current primer vector magnitude
%  pvdm = current primer derivative magnitude

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute state transition matrix at current time t

[~, ~, stm] = stm2(mu, t, ri, vi);

% evaluate primer vector fundamental equation

ppdot = stm * [pvi; pvdi];

% extract primer vector and primer derivative vector

pv = ppdot(1:3);

pvd = ppdot(4:6);

% compute primer vector magnitude

pvm = norm(pv);

% compute primer derivative magnitude

pvdm = dot(pv, pvd) / pvm;
