close all; clear; clc; addpath("..")

%% 2
R_earth = CelestialParameters.radius_earth;
mu_earth = CelestialParameters.gravityParameter_earth;

alt_Cygnss = 520; % km
r_Cygnss = R_earth + alt_Cygnss
v_Cygnss = sqrt(mu_earth / r_Cygnss)

alt_Gps = 20200; % km
r_Gps = R_earth + alt_Gps
v_Gps = sqrt(mu_earth / r_Gps)

v_GpsLos = v_Gps*r_cygnss/r_Gps
v_relMax = v_GpsLos + v_Cygnss