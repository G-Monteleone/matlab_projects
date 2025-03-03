function [uplvlh, pitch, yaw] = eci2lvlh (r, v, upeci)

% convert eci unit pointing vector to local
% vertical local horizontal coordinate system

% input

%  r     = eci position vector (kilometers)
%  v     = eci velocity vector (kilometers/second)
%  upeci = eci unit pointing vector (non-dimensional)

% output

%  uplvlh = lvlh unit pointing vector (non-dimensional)
%  pitch  = lvlh pitch angle (radians)
%           (-pi/2 <= pitch <= +pi/2)
%  yaw    = lvlh yaw angle (radians)
%           (0 <= yaw <= 2 pi)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tm = zeros(3, 3);

% angular momentum vector

hvm = cross(r, v);
   
hvms = dot(hvm, hvm);
   
hvmag = sqrt(hvms);

% unit angular momentum vector

hhat = hvm / hvmag;

rms = dot(r, r);
   
rmag = sqrt(rms);

% unit position vector

rhat = r / rmag;

hr = cross(hhat, rhat);

% compute components of eci-to-lvlh transformation matrix

for j = 1:1:3
    
    tm(1, j) = -hhat(j);
    
    tm(2, j) = hr(j);
    
    tm(3, j) = rhat(j);
    
end

% perform transformation

uplvlh = tm * upeci;

% compute lvlh angles (radians)

pitch = asin(uplvlh(3));

yaw = atan3(uplvlh(1), uplvlh(2));
