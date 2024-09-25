addpath('..')
clear
close all

HU = HypersonicsUtilities;

%% 1
densityTable = altitudeDensityTable';
calculatedAltitude = HU.airDensity(densityTable(1, :)*1000);
figure
plot(densityTable(2, :), densityTable(1, :), calculatedAltitude, densityTable(1, :))
legend("Tabulated", "Calculated")
title("Tabulated vs. Calculated Air Density vs. Altitude")
ylabel("Altitude (km)")
xlabel("Air Density (kg/m^3)")

figure
error = (calculatedAltitude - densityTable(2, :))./densityTable(2, :);
plot(error, densityTable(1, :))
title("Error Across Altitudes")
ylabel("Altitude (km)")
xlabel("Error ((calc - table)/table)")

%% 2
R_A = .05; % m
R_B = .9; % m
A = pi*R_B^2 % m^2
m = 3400; % kg
c_D = .4; % drag coefficient
k = 2e-4/sqrt(R_A)
v_0 = 4500; % m/s
h_0 = 100*1000; % m
gamma = deg2rad(-70); % flight path angle, radians

%% 2a
beta = HU.beta(m, gamma, c_D, A)
heating_max = HU.maxHeating(k, beta, v_0)
v_heatingMax = HU.velocityRatio_heatingMax*v_0
rho_heatingMax = HU.airDensity_maxHeating(beta)
h_heatingMax = HU.altitudeFromAirDensity(rho_heatingMax)

decel_max = HU.maxG(v_0, gamma)
rho_gMax = beta/2
h_gMax = HU.altitudeFromAirDensity(rho_gMax)
v_gMax = HU.velocityRatio_decelMax*v_0

%% 2b
% X = ballisticTrajectorySolver(1e-3, 30, v_0, gamma, h_0, m, c_D, A, false);
% XAltitude = [X(1:3, :); HU.altitude(X(4, :)); X(5, :)];
% Utilities.multiplot(XAltitude, ["velocity" "flight path angle" "altitude" "Earth angle"], ["time (s)" "velocity (m/s)" "\gamma (radians)" "altitude (m)" "\theta (radians)"])

% Y = [XAltitude(4, :); XAltitude(2:3, :); XAltitude(5, :)];
% Utilities.multiplotY(Y, ["velocity" "flight path angle" "Earth angle"], ["altitude (km)" "velocity (m/s)" "\gamma (radians)" "\theta (radians)"])

XSimplified = ballisticTrajectorySolver(1e-3, 28.5, v_0, gamma, h_0, m, c_D, A, true);