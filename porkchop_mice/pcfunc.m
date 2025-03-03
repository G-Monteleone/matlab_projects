function [dv1, dv2] = pcfunc (jdate1, jdate2)

% Lambert delta-v function

% required by porkchop.m

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu ip1 ip2 revmax

% time-of-flight (seconds)

taud = jdate2 - jdate1;

tof = taud * 86400.0;

% compute initial state vector

[ri, vi] = p2000(11, ip1, jdate1);

% compute final state vector

[rf, vf] = p2000(11, ip2, jdate2);
     
% solve Lambert's problem

sv1(1:3) = ri;

sv1(4:6) = vi;
    
sv2(1:3) = rf;

sv2(4:6) = vf;

[vito, vfto] = glambert(mu, sv1, sv2, tof, revmax);

% calculate departure delta-v vector

dv1(1) = vito(1) - vi(1);
dv1(2) = vito(2) - vi(2);
dv1(3) = vito(3) - vi(3);

% calculate arrival delta-v vector

dv2(1) = vf(1) - vfto(1);
dv2(2) = vf(2) - vfto(2);
dv2(3) = vf(3) - vfto(3);



