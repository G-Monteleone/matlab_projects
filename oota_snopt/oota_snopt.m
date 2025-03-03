% oota_snopt.m        February 20, 2024

% optimal impulsive orbital transfer analysis

% SNOPT version

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% global information

global dtr mu pvi pvdi

global top2eci sma ecc inc argper raan

global r1 v1 r2 v2 vt1 vt2 dv1 dv2

% angular conversion factors

dtr = pi / 180.0;

rtd = 180.0 / pi;

clc; home;

fprintf('\n program oota_snopt - SNOPT version\n');

fprintf('\n< optimal orbital transfer analysis >\n');

% request name of user-defined input data file

[filename, pathname] = uigetfile('*.in', 'please select an input data file');

% read simulation definition data file

[fid, mu, req, oev1, oev2, phi01, phi02, ...
    dphi1, dphi2, nphi1, nphi2] = read_oota(filename);

% load classical orbital elements of the initial orbit

sma(1) = oev1(1);
ecc(1) = oev1(2);
inc(1) = oev1(3) * dtr;
argper(1) = oev1(4) * dtr;
raan(1) = oev1(5) * dtr;

% load classical orbital elements of the final orbit

sma(2) = oev2(1);
ecc(2) = oev2(2);
inc(2) = oev2(3) * dtr;
argper(2) = oev2(4) * dtr;
raan(2) = oev2(5) * dtr;

% perform initialization

oota_iniz;

% convert initial search angles to radians

phi01 = phi01 * dtr;

phi02 = phi02 * dtr;

% convert initial search angle increments to radians

dphi1 = dphi1 * dtr;

dphi2 = dphi2 * dtr;

gdvmin = 1.0e99;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform two-dimensional grid search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of control variables

ncv = 2;

xg = zeros(ncv, 1);

xlwr = zeros(ncv, 1);

xupr = zeros(ncv, 1);

% number of mission constraints

nmc = 0;

flow = zeros(nmc, 1);
fupp = zeros(nmc, 1);

% bounds on objective function

flow(1) = 0.0e0;

fupp(1) = +Inf;

% solve the orbital TPBVP using SNOPT

xmul = zeros(ncv, 1);

xstate = zeros(ncv, 1);

fmul = zeros(nmc, 1);

fstate = zeros(nmc, 1);

tol = 1.0e-4;

for i = 1:1:nphi1

    for j = 1:1:nphi2

        % current initial guesses and bounds

        xg(1) = phi01 + dphi1 * (i - 1);

        xg(2) = phi02 + dphi2 * (j - 1);

        xlwr(1) = xg(1) - 2.0 * pi;

        xupr(2) = xg(2) + 2.0 * pi;

        % save current values

        xs1 = xg(1);

        xs2 = xg(2);

        dx = abs(xg(2) - xg(1));

        if ((abs(dx - pi) > tol) && (dx > tol) && (abs(dx - 2.0 * pi) > tol))

            % perform 2-dimensional minimization

            [xmin, ~, inform, xmul, fmul] = snopt(xg, xlwr, xupr, xmul, xstate, ...
                flow, fupp, fmul, fstate, 'oota_func1');

            % evaluate the current solution

            f = oota_func1(xmin);

            % check for global minimum

            if (f < gdvmin)

                % save current solution as the global minimum

                fs = f;

                gdvmin = fs;

                xsav(1) = xmin(1);

                xsav(2) = xmin(2);

            end

        end

    end

end

% compute characteristics of global optimum solution

f = oota_func1(xsav);

% compute eci delta-v vectors (kilometers/second)

dv1eci = top2eci * dv1';

dv2eci = top2eci * dv2';

% compute eci delta-v magnitudes (kilometers/second)

dv1ecim = norm(dv1eci);

dv2ecim = norm(dv2eci);

% ------------------------------------------------------------
% compute eci state vectors (kilometers and kilometers/second)
% ------------------------------------------------------------

% initial orbit - prior to first impulse

r1eci = top2eci * r1';

v1eci = top2eci * v1';

oev1 = eci2orb1 (mu, r1eci, v1eci);

upeci = dv1eci / dv1ecim;

[~, pitch1, yaw1] = eci2lvlh (r1eci, v1eci, upeci);

fprintf('\ninitial orbit - prior to the first impulse');
fprintf('\n------------------------------------------\n');

oeprint1(mu, oev1, 1);

[r, v] = orb2eci(mu, oev1);

svprint(r, v);

% transfer orbit - after the first impulse

rt1eci = top2eci * r1';

vt1eci = top2eci * vt1';

oevt1 = eci2orb1 (mu, rt1eci, vt1eci);

fprintf('\ntransfer orbit - after the first impulse');
fprintf('\n----------------------------------------\n');

oeprint1(mu, oevt1, 1);

[r, v] = orb2eci(mu, oevt1);

svprint(r, v);

% transfer orbit - prior to second impulse

rt2eci = top2eci * r2';

vt2eci = top2eci * vt2';

oevt2 = eci2orb1 (mu, rt2eci, vt2eci);

upeci = dv2eci / dv2ecim;

[uplvlh, pitch2, yaw2] = eci2lvlh (rt2eci, vt2eci, upeci);

fprintf('\ntransfer orbit - prior to the second impulse');
fprintf('\n--------------------------------------------\n');

oeprint1(mu, oevt2, 1);

[r, v] = orb2eci(mu, oevt2);

svprint(r, v);

% final orbit - after the second impulse

r2eci = top2eci * r2';

v2eci = top2eci * v2';

oev2 = eci2orb1 (mu, r2eci, v2eci);

fprintf('\nfinal orbit - after the second impulse');
fprintf('\n--------------------------------------\n');

oeprint1(mu, oev2, 1);

[r, v] = orb2eci(mu, oev2);

svprint(r, v);

fprintf('\nECI delta-v vectors, magnitudes and LVLH angles');
fprintf('\n-----------------------------------------------\n');

fprintf('\ndelta-v1x         %12.4f meters/second', 1000 * dv1eci(1));
fprintf('\ndelta-v1y         %12.4f meters/second', 1000 * dv1eci(2));
fprintf('\ndelta-v1z         %12.4f meters/second', 1000 * dv1eci(3));
fprintf('\n\ndelta-v1          %12.4f meters/second\n', 1000 * dv1ecim);

fprintf('\nLVLH pitch angle  %12.4f degrees', rtd * pitch1);
fprintf('\nLVLH yaw angle    %12.4f degrees\n', rtd * yaw1);

fprintf('\ndelta-v2x         %12.4f meters/second', 1000 * dv2eci(1));
fprintf('\ndelta-v2y         %12.4f meters/second', 1000 * dv2eci(2));
fprintf('\ndelta-v2z         %12.4f meters/second', 1000 * dv2eci(3));
fprintf('\n\ndelta-v2          %12.4f meters/second\n', 1000 * dv2ecim);

fprintf('\nLVLH pitch angle  %12.4f degrees', rtd * pitch2);

fprintf('\nLVLH yaw angle    %12.4f degrees\n', rtd * yaw2);

fprintf('\ntotal delta-v     %12.4f meters/second\n', ...
    1000 * (dv1ecim + dv2ecim));

% compute flight time on transfer orbit (seconds)

dtof = tof1(mu, oevt1(1), oevt1(2), oevt1(6), oevt2(6));

fprintf('\ntransfer time    %12.4f seconds', dtof);

fprintf('\n                 %12.4f minutes', dtof / 60.0);

fprintf('\n                 %12.4f hours\n', dtof / 3600.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display trajectory graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of graphic data points

npts = 300;

% orbital periods of initial and final orbits (seconds)

period1 = 2.0 * pi * oev1(1) * sqrt(oev1(1) / mu);

period3 = 2.0 * pi * oev2(1) * sqrt(oev2(1) / mu);

deltat1 = period1 / npts;

simtime1 = -deltat1;

deltat3 = period3 / npts;

simtime3 = -deltat3;

% allocate graphic data points

rp1_x = zeros(npts + 1, 1);

rp1_y = zeros(npts + 1, 1);

rp1_z = zeros(npts + 1, 1);

rp2_x = zeros(npts + 1, 1);

rp2_y = zeros(npts + 1, 1);

rp2_z = zeros(npts + 1, 1);

rp3_x = zeros(npts + 1, 1);

rp3_y = zeros(npts + 1, 1);

rp3_z = zeros(npts + 1, 1);

% compute graphic data points for initial and final orbits

for i = 1:1:npts + 1

    simtime1 = simtime1 + deltat1;

    simtime3 = simtime3 + deltat3;

    % compute initial orbit "normalized" position vector

    [rwrk, ~] = twobody2(mu, simtime1, r1eci, v1eci);

    rp1_x(i) = rwrk(1) / req;

    rp1_y(i) = rwrk(2) / req;

    rp1_z(i) = rwrk(3) / req;

    % compute final orbit "normalized" position vector

    [rwrk, vwrk] = twobody2(mu, simtime3, r2eci, v2eci);

    rp3_x(i) = rwrk(1) / req;

    rp3_y(i) = rwrk(2) / req;

    rp3_z(i) = rwrk(3) / req;

end

% compute graphic data points for transfer orbit

deltat2 = dtof / npts;

simtime2 = -deltat2;

for i = 1:1:npts + 1

    simtime2 = simtime2 + deltat2;

    % compute initial orbit "normalized" position vector

    [rwrk, ~] = twobody2(mu, simtime2, rt1eci, vt1eci);

    rp2_x(i) = rwrk(1) / req;

    rp2_y(i) = rwrk(2) / req;

    rp2_z(i) = rwrk(3) / req;

end

figure(1);

% create axes vectors

xaxisx = [1 1.5];
xaxisy = [0 0];
xaxisz = [0 0];

yaxisx = [0 0];
yaxisy = [1 1.5];
yaxisz = [0 0];

zaxisx = [0 0];
zaxisy = [0 0];
zaxisz = [1 1.5];

figure (1);

hold on;

grid on;

% plot spherical earth

[x, y, z] = sphere(24);

h = surf(x, y, z);

colormap([127/255 1 222/255]);

set (h, 'edgecolor', [1 1 1]);

% plot coordinate system axes

plot3(xaxisx, xaxisy, xaxisz, '-g', 'LineWidth', 1);

plot3(yaxisx, yaxisy, yaxisz, '-r', 'LineWidth', 1);

plot3(zaxisx, zaxisy, zaxisz, '-b', 'LineWidth', 1);

% plot initial orbit

plot3(rp1_x, rp1_y, rp1_z, '-r', 'LineWidth', 1.5);

% plot final orbit

plot3(rp3_x, rp3_y, rp3_z, '-g', 'LineWidth', 1.5);

% plot transfer orbit

plot3(rp2_x, rp2_y, rp2_z, '-b', 'LineWidth', 2.0);

% label eci axes in earth radii

xlabel('X coordinate (ER)', 'FontSize', 12);

ylabel('Y coordinate (ER)', 'FontSize', 12);

zlabel('Z coordinate (ER)', 'FontSize', 12);

% label orbital locations of initial and final impulses

plot3(rp1_x(1), rp1_y(2), rp1_z(3), 'ob', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'b');

plot3(rp2_x(end), rp2_y(end), rp2_z(end), 'ob', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'b');

title('Optimal Orbit Transfer Analysis', ...
    'Initial, Transfer and Final Orbits', 'FontSize', 16);

axis equal;

view(50, 20);

rotate3d on;

%%%%%%%%%%%%%%%%%%%%%%%%
% primer vector analysis
%%%%%%%%%%%%%%%%%%%%%%%%

if (dv1ecim > 0.0 && dv2ecim > 0.0)
    
    % perform primer vector initialization
    
    pviniz(dtof, rt1eci', vt1eci', dv1eci', dv2eci');
    
    % number of graphic data points
    
    npts = 300;
    
    x1 = zeros(npts + 1, 1);
        
    y1 = zeros(npts + 1, 1);
        
    y2 = zeros(npts + 1, 1);

    % plot behavior of primer vector magnitude
    
    dt = dtof / npts;
    
    for i = 1:1:npts + 1
        
        t = (i - 1) * dt;
        
        if (t == 0)
            
            % initial value of primer magnitude and derivative
            
            pvm = norm(pvi);
            
            pvdm = dot(pvi, pvdi) / pvm;
            
        else
            
            % primer vector and derivative magnitudes at time t
            
            [pvm, pvdm] = pvector(rt1eci, vt1eci, t);
            
        end
        
        % load data array

        x1(i) = t / 60.0;
        
        y1(i) = pvm;
        
        y2(i) = pvdm;
        
    end
    
    % create primer graphics
    
    figure(2);
    
    hold on;
    
    plot(x1, y1, '-r');
    
    plot(x1, y1, '.r');
    
    title('Optimal Orbit Transfer Analysis', ...
        'Primer Vector Analysis', 'FontSize', 16);
    
    xlabel('simulation time (minutes)', 'FontSize', 12);
    
    ylabel('primer vector magnitude', 'FontSize', 12);
    
    grid;
    
    % plot behavior of magnitude of primer derivative
    
    figure(3);
    
    hold on;
    
    plot(x1, y2, '-r');
    
    plot(x1, y2, '.r');
    
    title('Optimal Orbit Transfer Analysis', ...
        'Primer Vector Analysis', 'FontSize', 16);
    
    xlabel('simulation time (seconds)', 'FontSize', 12);
    
    ylabel('primer derivative magnitude', 'FontSize', 12);
    
    grid;
    
end