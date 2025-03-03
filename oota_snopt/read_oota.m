function [fid, mu, req, oev1, oev2, phi01, phi02, ...
    dphi1, dphi2, nphi1, nphi2] = read_oota(filename)

% read orbital elements and simulation control data file

% input

%  filename = name of oota data file

% output

%  fid = file id

%  oev1 = initial orbit orbital elements array
%  oev2 = final orbit orbital elements array

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open data file

fid = fopen(filename, 'r');

% check for file open error

if (fid == -1)
    
    fprintf('\n\n  error: cannot find this file!!');
    
    pause
    
    return;
    
end

% read 81 lines of data

for i = 1:1:81
    
    cline = fgetl(fid);
    
    switch i
        
        case 8
            
            % gravitational constant (kilometers^3/seconds^2)

            mu = str2double(cline);
        
        case 11
            
            % earth equatorial radius (kilometers)

            req = str2double(cline);
            
        case 19
            
            % semimajor axis of the initial orbit (kilometers)

            oev1(1) = str2double(cline);
            
        case 23
            
            % orbital eccentricity of the initial orbit

            oev1(2) = str2double(cline);
            
        case 27
            
            % orbital inclination of the initial orbit (degrees)

            oev1(3) = str2double(cline);
            
        case 31
            
            % argument of perigee of the initial orbit (degrees)

            oev1(4) = str2double(cline);
            
        case 35
            
            % raan of the initial orbit (degrees)

            oev1(5) = str2double(cline);
            
        case 43
            
            % semimajor axis of the final orbit (kilometers)

            oev2(1) = str2double(cline);
            
        case 47
            
            % orbital eccentricity of the final orbit

            oev2(2) = str2double(cline);
            
        case 51
            
            % orbital inclination of the final orbit (degrees)

            oev2(3) = str2double(cline);
            
        case 55
            
            % argument of perigee of final orbit (degrees)

            oev2(4) = str2double(cline);
            
        case 59
            
            % raan of final orbit (degrees)

            oev2(5) = str2double(cline);
            
        case 66
            
            % initial orbit true anomaly at which to begin search (degrees)

            phi01 = str2double(cline);
            
        case 69
            
            % final orbit true anomaly at which to begin search (degrees)

            phi02 = str2double(cline);
            
        case 72
            
            % initial orbit true anomaly search increment (degrees)

            dphi1 = str2double(cline);
            
        case 75
            
            % final orbit true anomaly search increment (degrees)

            dphi2 = str2double(cline);
            
        case 78
            
            % number of initial orbit true anomaly search intervals

            nphi1 = str2double(cline);
            
        case 81
            
            % number of final orbit true anomaly search intervals
            
            nphi2 = str2double(cline);
            
    end
end

fclose(fid);

