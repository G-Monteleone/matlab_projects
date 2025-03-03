clc
clear all
format long g

re = 6378.1363; % [km]
r1 = [ 2.500000,    0.000000 ,   0.000000]*re; % [km]
r2 = [ 1.9151111,   1.6069690,   0.000000]*re; % [km]
tof = 76*60; % [s]

[V1, V2] = LAMBERTBATTIN(r1, r2, 'pro', tof)

