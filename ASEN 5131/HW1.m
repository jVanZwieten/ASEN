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
h_0 = 100; % km
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
w = m*HU.earthSurfaceGravityAcceleration;

dt = .001; % s
t_f = 28.5; % s
steps = t_f/dt + 1; % +1 accounts for t = 0;
X = [zeros(5, length(steps))]; % [t, v, gamma, r, theta]
X(:, 1) = [0; v_0; gamma; h_0*1000 + HU.earthRadius; 0];
for i = 2:steps
    X(:, i) = TimeStep(X(:, i - 1), dt, w, c_D, A);
end

XAltitude = [X(1:3, :); HU.altitude(X(4, :)); X(5, :)];
Utilities.multiplot(XAltitude, ["velocity" "flight path angle" "altitude" "Earth angle"], ["time (s)" "velocity (m/s)" "\gamma (radians)" "altitude (m)" "\theta (radians)"])

Y = [XAltitude(4, :); XAltitude(2:3, :); XAltitude(5, :)];
Utilities.multiplotY(Y/1000, ["velocity" "flight path angle" "Earth angle"], ["altitude (km)" "velocity (m/s)" "\gamma (radians)" "\theta (radians)"])

function X_next = TimeStep(X, dt, weight, dragCoefficient, area)
    h = HypersonicsUtilities.altitude(X(4));
    rho = HypersonicsUtilities.airDensity(h);
    d = HypersonicsUtilities.drag(dragCoefficient, area, rho, X(2));

    X_next = zeros(5, 1);
    X_next(1) = X(1) + dt;
    X_next(2) = X(2) + dt*HypersonicsUtilities.acceleration(0, weight, d, X(3));
    X_next(3) = X(3) + dt*HypersonicsUtilities.flightPathAngleRate(X(2), 0, weight, X(4), X(3));
    X_next(4) = X(4) + dt*HypersonicsUtilities.radialDistanceRate(X(2), X(3));
    X_next(5) = X(5) + dt*HypersonicsUtilities.earthAngleRate(X(2), X(3), X(4));
end