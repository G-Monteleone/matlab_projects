function [f, g] = iptofunc (x)

% delta-v objective function and mission constraints

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global smu ip1 ip2 jdtdb0 eq2000

global otype imcon dv_depart dv_arrive

global jdtdb1 jdtdb2 ri vi rf vf vito vfto sv2

dv_depart = zeros(3, 1);

dv_arrive = zeros(3, 1);

% current tdb julian dates

jdtdb1 = x(1) + jdtdb0;

jdtdb2 = x(2) + jdtdb0;

% time-of-flight (seconds)

tof = 86400.0 * (jdtdb2 - jdtdb1);

% compute initial state vector

[ri, vi] = p2000_ecl(ip1, jdtdb1);

% compute final state vector

[rf, vf] = p2000_ecl(ip2, jdtdb2);

% solve Lambert's problem

sv1(1:3) = ri;

sv1(4:6) = vi;

sv2(1:3) = rf;

sv2(4:6) = vf;

[vito, vfto] = glambert(smu, sv1, sv2, tof, 0.0);

% calculate departure delta-v (kilometers/second)

dv_depart(1) = vito(1) - vi(1);
dv_depart(2) = vito(2) - vi(2);
dv_depart(3) = vito(3) - vi(3);

% calculate arrival delta-v (kilometers/second)

dv_arrive(1) = vf(1) - vfto(1);
dv_arrive(2) = vf(2) - vfto(2);
dv_arrive(3) = vf(3) - vfto(3);

if (imcon == 1)
    
    f = zeros(5, 1);
    
    % -----------------------------------
    % compute current mission constraints
    % -----------------------------------
    
    % C3L (kilometers^2/second^2)
    
    f(2) = norm(dv_depart) * norm(dv_depart);
    
    % DLA (radians)
    
    if (ip1 == 3)
        
        % Earth is departure planet
        
        dv_eq = eq2000 * dv_depart;
        
        f(3) = 0.5 * pi - acos(dv_eq(3) / norm(dv_eq));
        
    else
        
        f(3) = 0.5 * pi - acos(dv_depart(3) / norm(dv_depart));
        
    end
    
    % TOF (days)
    
    f(4) = jdtdb2 - jdtdb1;
    
    % arrival v-infinity (kilometers/second)
    
    f(5) = norm(dv_arrive);
    
else
    
    f = zeros(1, 1);
    
end

% load scalar objective function

switch otype
    
    case 1
        
        % launch
        
        f(1) = norm(dv_depart);
        
    case 2
        
        % arrival
        
        f(1) = norm(dv_arrive);
        
    case 3
        
        % launch + arrival
        
        f(1) = norm(dv_depart) + norm(dv_arrive);
        
end

% no derivatives

g = [];
