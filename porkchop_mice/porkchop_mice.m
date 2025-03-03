% porkchop_mice.m        January 16, 2022

% MATLAB script for creating ballistic Earth-to-Mars
% interplanetary pork chop plots

% DE421 ephemeris and MICE routines

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global mu iephem ephname km ip1 ip2 revmax

rtd = 180.0 / pi;

% gravitational constant of the sun (km^3/sec^2)

mu = 132712441933.0;

% initialize de421 ephemeris

ephname = 'de421.bsp';

iephem = 1;

km = 1;

row = 0;
col = 0;

clc; home;

fprintf('\nprogram porkchop_mice - interplanetary pork chop plots');

fprintf('\n\nnominal launch date\n');

[month, day, year] = getdate;

JD_dep = julian(month, day, year);

fprintf('\nnominal arrival date\n');

[month, day, year] = getdate;

JD_arr = julian(month, day, year);

while(1)
    
   fprintf('\nplease input the launch date span in days\n');
   
   spanx = input('? ');
   
   if (spanx > 0)
   
      break;
      
   end
   
end 

while(1)
    
   fprintf('\nplease input the arrival date span in days\n');
   
   spany = input('? ');
   
   if (spany > 0)
   
      break;
      
   end
   
end

while(1)
    
   fprintf('\nplease input the step size in days\n');
   
   step = input('? ');
   
   if (step > 0)
   
      break;
      
   end
   
end

fprintf('\nplease input the launch energy contour levels in km^2/sec^2 ([] for defaults)\n');
    
c3_levels = input('? ');

if (size(c3_levels, 2) == 0)
    
   c3_levels = [6,7,8,9,10,11,12,13,14,15,16,20,25,30,35,40,45,50];
   
end
   
fprintf('\nplease input the arrival v-infinity contour levels in km/sec ([] for defaults)\n');
    
vinf_levels = input('? ');

if (size(vinf_levels, 2) == 0)
    
   vinf_levels = [1.5,1.8,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.5,4.0,5.0,6.0,7.0,8.0];
   
end

fprintf('\nplease input the launch and arrival declination contour levels in degrees ([] for defaults)\n');
    
dla_levels = input('? ');

if (size(dla_levels, 2) == 0)
    
   dla_levels = [-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30];
   
end

fprintf('\nplease input the launch and arrival right ascension contour levels in degrees ([] for defaults)\n');
    
rla_levels = input('? ');

if (size(rla_levels, 2) == 0)
    
   rla_levels = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330];
   
end

fprintf('\nplease input the time-of-flight contour levels in days ([] for defaults)\n');
    
tof_levels = input('? ');

if (size(tof_levels, 2) == 0)
    
   tof_levels = [100, 150, 200, 250, 300, 350, 400, 450];
   
end

fprintf('\nplease input the total delta-v contour levels in kilometers/second ([] for defaults)\n');
    
dvt_levels = input('? ');

if (size(dvt_levels, 2) == 0)
    
   dvt_levels = [5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12.0, 13.0, 14.0, 15.0];
   
end

% launch planet is Earth

ip1 = 3;

% arrival planet is Mars

ip2 = 4;

% restrict to type I and II trajectories

revmax = 0;

% compute C3L, v-infinity, DLA and RLA data

for i = JD_dep - spanx: step: JD_dep + spanx
    
    row = row + 1;
    
    JDi = i;
    
    col = 0;
    
    for j = JD_arr - spany: step: JD_arr + spany
        
        col = col + 1;
        
        JDf = j;
        
        con_x(row, col) = JDi - JD_dep;
        
        con_y(row, col) = JDf - JD_arr;
                
        [dv1, dv2] = pcfunc (JDi, JDf);
                      
        c3l = norm(dv1) * norm(dv1);
        
        c3_dep(row, col) = c3l;
               
        vinf_arr(row, col) = norm(dv2);
        
        if (c3l <= 50)
            
           % compute orientation of the departure hyperbola
  
           decl_dep(row, col) = 90.0 - rtd * acos(dv1(3) / norm(dv1));

           rasc_dep(row, col) = rtd * atan3(dv1(2), dv1(1));
           
        else
            
           decl_dep(row, col) = NaN;

           rasc_dep(row, col) = NaN;
           
        end
       
        % compute arrival dla and rla
        
        tmatrix = mme2000 (JD_dep + i);
        
        dvmmee = tmatrix * (-dv2');
        
        decl_arr(row, col) = 90.0 - rtd * acos(dvmmee(3) / norm(dv2));

        rasc_arr(row, col) = rtd * atan3(dvmmee(2), dvmmee(1));
        
        % compute flight time in days
        
        tof(row, col) = JDf - JDi;
        
        % compute total delta-v in kilometers/second
        
        dvtotal(row, col) = norm(dv1) + norm(dv2);
        
    end
    
end

% create string representations of launch and arrival dates

[cdstr1, utstr1] = jd2str(JD_dep);

[cdstr2, utstr2] = jd2str(JD_arr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3l and arrival v-infinity plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);

[c1, h1] = contour(con_x, con_y, c3_dep, c3_levels, 'r');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, vinf_arr, vinf_levels, 'b');
 
clabel(c2, h2);

title({'Ballistic Earth-to-Mars Trajectory'; 'launch C3 and arrival v-infinity'}, 'FontSize', 16');

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Mars arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'C3L (kilometers^2/seconds^2)', 'arrival v-infinity (kilometers/second)', 'Location', 'northwest');

print -depsc -tiff -r300 c3l_vinf.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3L and DLA plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);

[c1, h1] = contour(con_x, con_y, decl_dep, dla_levels, 'b');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, c3_dep, c3_levels, 'r');
 
clabel(c2, h2);

title({'Ballistic Earth-to-Mars Trajectory'; 'launch C3 and DLA'}, 'FontSize', 16');

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Mars arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'DLA (degrees)', 'C3L (kilometers^2/second^2)', 'Location', 'northwest');

print -depsc -tiff -r300 c3l_dla.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3L and RLA plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);

[c1, h1] = contour(con_x, con_y, rasc_dep, rla_levels, 'b');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, c3_dep, c3_levels, 'r');
 
clabel(c2, h2);

title({'Ballistic Earth-to-Mars Trajectory'; 'launch C3 and RLA'}, 'FontSize', 16');

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Mars arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'RLA (degrees)', 'C3L (kilometers^2/second^2)', 'Location', 'northwest');

print -depsc -tiff -r300 c3l_rla.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined C3L and flight time plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);

[c1, h1] = contour(con_x, con_y, tof, tof_levels, 'b');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, c3_dep, c3_levels, 'r');
 
clabel(c2, h2);

title({'Ballistic Earth-to-Mars Trajectory'; 'launch C3 and flight time'}, 'FontSize', 16);

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Mars arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'TOF (days)', 'C3L (kilometers^2/second^2)', 'Location', 'northwest');

print -depsc -tiff -r300 c3l_tof.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create combined arrival v-infinity and dla plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5);

[c1, h1] = contour(con_x, con_y, vinf_arr, vinf_levels, 'b');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, decl_arr, dla_levels, 'r');
 
clabel(c2, h2);

title({'Ballistic Earth-to-Mars Trajectory'; 'arrival v-infinity and DLA'}, 'FontSize', 16);

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Mars arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'v-infinity (kilometers/second)', 'DLA (degrees)', 'Location', 'northwest');

print -depsc -tiff -r300 vinf_dla.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create arrival v-infinity and rla plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6);

[c1, h1] = contour(con_x, con_y, vinf_arr, vinf_levels, 'b');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, rasc_arr, rla_levels, 'r');
 
clabel(c2, h2);

title({'Ballistic Earth-to-Mars Trajectory'; 'arrival v-infinity and RLA'}, 'FontSize', 16);

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Mars arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'v-infinity (kilometers/second)', 'RLA (degrees)', 'Location', 'northwest');

print -depsc -tiff -r300 vinf_rla.eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create toal delta-v and flight time plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7);

[c1, h1] = contour(con_x, con_y, tof, tof_levels, 'b');

clabel(c1, h1);

hold on;

[c2, h2] = contour(con_x, con_y, dvtotal, dvt_levels, 'r');
 
clabel(c2, h2);

title({'Ballistic Earth-to-Mars Trajectory'; 'total delta-v and flight time'}, 'FontSize', 16);

xlabel(['days relative to Earth departure date of ', cdstr1]);

ylabel(['days relative to Mars arrival date of ', cdstr2]);

legend([h1(1), h2(1)], 'TOF (days)', 'DVT (kilometers/second)', 'Location', 'northwest');

print -depsc -tiff -r300 dvt_tof.eps;

% unload ephemeris

cspice_unload('de421.bsp')
