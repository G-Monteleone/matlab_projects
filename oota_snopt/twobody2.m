function [rf, vf] = twobody2 (mu, tau, ri, vi)

% solve the two body initial value problem

% Shepperd's method

% input

%  mu  = gravitational constant (kilometers^3/seconds^2)
%  tau = propagation time interval (seconds)
%  ri  = initial eci position vector (kilometers)
%  vi  = initial eci velocity vector (kilometers/second)

% output

%  rf = final eci position vector (kilometers)
%  vf = final eci velocity vector (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tolerance = 1.0e-12;

u = 0.0;

% maximum number of iterations

imax = 20;

umax = +realmax;

umin = -realmax;

orbits = 0;

tdesired = tau;

threshold = tolerance * abs(tdesired);

r0 = norm(ri);

n0 = dot(ri, vi);

beta = 2.0 * (mu / r0) - dot(vi, vi);

if (beta ~= 0.0)

    umax = +1.0 / sqrt(abs(beta));
    
    umin = -1.0 / sqrt(abs(beta));

end

if (beta > 0)

    orbits = beta * tau - 2.0 * n0;
    
    orbits = 1.0 + (orbits * sqrt(beta)) / (pi * mu);
    
    orbits = floor (orbits / 2);
end

for i = 1:1:imax

    q = beta * u * u;
    
    q = q / (1.0 + q);

    n =  0;
    
    r = 1.0;
    
    l = 1.0;
    
    s = 1.0;
    
    d = 3.0;
    
    gcf = 1.0;
    
    k = -5;
    
    gold = 0.0;

    while (gcf ~= gold)

        k = -k;
        
        l = l + 2.0;
        
        d = d + 4.0 * l;
        
        n = n + (1.0 + k) * l;
        
        r = d / (d - n * r * q);
        
        s = (r - 1.0) * s;

        gold = gcf;
        
        gcf  = gold + s;

    end

    h0 = 1.0 - 2.0 * q;
    
    h1 = 2.0 * u * (1.0 - q);

    u0 = 2.0 * h0 * h0 - 1.0;
    
    u1 = 2.0 * h0 * h1;
    
    u2 = 2.0 * h1 * h1;
    
    u3 = 2.0 * h1 * u2 * gcf / 3.0;

    if (orbits ~= 0)

        u3 = u3 + 2.0 * pi * orbits / (beta * sqrt(beta));

    end

    r1 = r0 * u0 + n0 * u1 + mu * u2;
    
    dt = r0 * u1 + n0 * u2 + mu * u3;
    
    slope = 4.0 * r1 / (1.0 + beta * u * u);
    
    terror = tdesired - dt;

    if (abs (terror) < threshold)

        break;

    end
    
    if ((i > 1) && (u  == uold))

        break;

    end
    
    if ((i > 1) && (dt == dtold))

        break;

    end

    uold  = u;
    
    dtold = dt;
    
    ustep = terror / slope;

    if (ustep > 0.0)

        umin = u;
        
        u = u + ustep;
        
        if (u > umax)

            u = (umin + umax) / 2.0;

        end

    else

        umax = u;
        
        u = u + ustep;
        
        if (u < umin)

            u = (umin + umax) / 2.0;

        end
    end

    if (i == imax)

        fprintf('\n\nmax iterations in twobody2 function');

        pause
    end

end

f = 1.0 - (mu / r0) * u2;

gg = 1.0 - (mu / r1) * u2;

g  =  r0 * u1 + n0 * u2;

ff = -mu * u1 / (r0 * r1);

% final position and velocity vectors (kilometers and kilometers/second)

rf = f  * ri + g  * vi;

vf = ff * ri + gg * vi;





