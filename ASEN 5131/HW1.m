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
dt = 1e-3; % s
t_f = 28.5; % s

X = ballisticTrajectorySolver(dt, t_f, v_0, gamma, h_0, m, c_D, A, false);
plotBallisticSim(X, 'Ballistic Simulation')
exportBallisticSim(X, 'ballisticSim.dat')

XSimplified = ballisticTrajectorySolver(dt, t_f, v_0, gamma, h_0, m, c_D, A, true);
plotBallisticSim(XSimplified, 'Ballistic Simulation Simplified')
exportBallisticSim(XSimplified, 'ballisticSimSimplified.dat')

%% 2c
[maxQ_Simplified h_maxQSimplified v_maxQSimplified] = maxQ(XSimplified, k)
[maxG_Simplified h_maxGSimplified v_maxGSimplified] = maxG(XSimplified, c_D, A, m)

%% 2d
[maxQ_Sim h_maxQSim v_maxQSim] = maxQ(X, k)
[maxG_Sim h_maxGSim v_maxGSim] = maxG(X, c_D, A, m)

%% 2e
plotVelocityTime(X, 'Ballistic Simulation, Unsimplified')
plotVelocityTime(XSimplified, 'Ballistic Simulation Simplified')

function plotBallisticSim(X, title)
    XAltitude = [X(1:3, :); HypersonicsUtilities.altitude(X(4, :)); X(5, :)];
    Utilities.multiplot(XAltitude, ["velocity" "flight path angle" "altitude" "Earth angle"], ["time (s)" "velocity (m/s)" "\gamma (radians)" "altitude (m)" "\theta (radians)"], title + " time plots")

    Y = [XAltitude(4, :); XAltitude(2:3, :); XAltitude(5, :)];
    Utilities.multiplotY(Y, ["velocity" "flight path angle" "Earth angle"], ["altitude (km)" "velocity (m/s)" "\gamma (radians)" "\theta (radians)"], title + " altitude profiles")
end

function exportBallisticSim(X, filename)
    format = '%12.5g';

    time = sprintfc('%6.4f', X(1, :)');
    velocity = sprintfc(format, X(2, :)');
    flightPathAngle = sprintfc(format, X(3, :)');
    altitude = sprintfc(format, HypersonicsUtilities.altitude(X(4, :))');
    earthAngle = sprintfc(format, X(5, :)');
    
    t = table(time, velocity, flightPathAngle, altitude, earthAngle);
    writetable(t, filename)
end

function [maxQ, h_maxQ, v_maxQ] = maxQ(X, k)
    Q = HypersonicsUtilities.heatTransfer(k, HypersonicsUtilities.airDensity(HypersonicsUtilities.altitude(X(4, :))), X(2, :));
    [maxQ i] = max(Q);
    h_maxQ = HypersonicsUtilities.altitude(X(4, i));
    v_maxQ = X(2, i);
end

function [maxG, h_maxG, v_maxG] = maxG(X, c_D, A, m)
    G = HypersonicsUtilities.acceleration(0, HypersonicsUtilities.drag(c_D, A, HypersonicsUtilities.airDensity(HypersonicsUtilities.altitude(X(4, :))), X(2, :)), m, X(3, :))/HypersonicsUtilities.earthSurfaceGravityAcceleration;
    [maxG i] = min(G);
    h_maxG = HypersonicsUtilities.altitude(X(4, i));
    v_maxG = X(2, i);
end

function plotVelocityTime(X, figureTitle)
    figure
    t = tiledlayout(2, 1);
    title(t, figureTitle)
    altitude = HypersonicsUtilities.altitude(X(4, :))/1000;

    nexttile
    plot(X(2, :), altitude)
    ylabel('altitude (km)')
    xlabel('velocity (m/s)')
    title('velocity profile')

    nexttile
    plot(X(1, :), altitude)
    ylabel('altitude (km)')
    xlabel('time (s)')
    title('time profile')
end