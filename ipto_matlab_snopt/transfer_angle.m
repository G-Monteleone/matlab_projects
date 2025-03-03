function tangle = transfer_angle(r1, r2)

% calculate the transfer or "included" angle between
% two coplanar position vectors

% input

%  r1 = first position vector
%  r2 = second position vector

% output

%  tangle = transfer angle (radians)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unit "normal" vector

xkhat(1) = 0.0;
xkhat(2) = 0.0;
xkhat(3) = 1.0;

% compute unit position vectors

rhat1 = r1 / norm(r1);

rhat2 = r2 / norm(r2);

% cosine of transfer angle

ctangle = dot(rhat1, rhat2);

r1xr2 = cross(r1, r2);

vdotwrk = dot(r1xr2, xkhat);

% sine of transfer angle

stangle = sign(vdotwrk) * sqrt(1.0 - ctangle^2);

% transfer angle (radians)

tangle = atan3(stangle, ctangle);

end

