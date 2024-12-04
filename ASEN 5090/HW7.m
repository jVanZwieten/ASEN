close all; clear; clc; addpath("..")

%% 2
R_earth = CelestialParameters.radius_earth/1000; % km
mu_earth = CelestialParameters.gravityParameter_earth;

alt_Cygnss = 520; % km
r_Cygnss = R_earth + alt_Cygnss
v_Cygnss = sqrt(mu_earth / r_Cygnss)

alt_Gps = 20200; % km
r_Gps = R_earth + alt_Gps
v_Gps = sqrt(mu_earth / r_Gps)

v_GpsLos = v_Gps*r_Cygnss/r_Gps
v_relMax = v_GpsLos + v_Cygnss

c = Utilities.speedOfLight/1000; % km/s
f_L1 = 1575.42e6; % Hz
f_doppler = f_L1*v_relMax/c