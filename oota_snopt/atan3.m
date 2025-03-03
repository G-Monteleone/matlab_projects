function y = atan3 (a, b)

% four quadrant inverse tangent

% input

%  a = sine of angle
%  b = cosine of angle

% output

%  y = angle (radians; 0 =< c <= 2 * pi)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 1.0e-10;

pidiv2 = 0.5 * pi;

if (abs(a) < epsilon)

   y = (1.0 - sign(b)) * pidiv2;

   return;

else

   c = (2.0 - sign(a)) * pidiv2;

end

if (abs(b) < epsilon)

   y = c;

   return;

else

   y = c + sign(a) * sign(b) * (abs(atan(a / b)) - pidiv2);

end
