clc; clear; close all;
format long e
addpath('..')

Gtil = CelestialParameters.universalGravitationalConstant;

%% 1a
GMtil_earth = CelestialParameters.gravityParameter_earth;
Mtil_earth = GMtil_earth/Gtil

GMtil_luna = CelestialParameters.gravityParameter_luna;
Mtil_luna = GMtil_luna/Gtil

mStar_earthLuna = CR3BPUtilities.characteristicMass(Mtil_earth, Mtil_luna)
mu_earthLuna = CR3BPUtilities.massRatio(Mtil_earth, Mtil_luna)

lStar_earthLuna = CelestialParameters.semiMajorAxis_luna;

tStar_earthLuna = CR3BPUtilities.characteristicTime(lStar_earthLuna, mStar_earthLuna)

GMtil_sol = CelestialParameters.gravityParameter_sol;
Mtil_sol = GMtil_sol/Gtil

mStar_solEarth = CR3BPUtilities.characteristicMass(Mtil_sol, Mtil_earth)
mu_solEarth = CR3BPUtilities.massRatio(Mtil_sol, Mtil_earth)

lStar_solEarth = CelestialParameters.semiMajorAxis_earth;

tStar_solEarth = CR3BPUtilities.characteristicTime(lStar_solEarth, mStar_solEarth)

%% 2bi
xVec_0 = [0.98; 0; 0; 0; 1.2; 0];
tSpan = [0 2];
tolerancePower = -12;
[t_i, X_i] = integrateCr3bp(xVec_0, mu_earthLuna, tSpan, tolerancePower);

figure; plot3(X_i(:, 1), X_i(:, 2), X_i(:, 3), 'DisplayName', 'Trajectory');
plotCr3bp(X_i, mu_earthLuna);

%% 2bii
xVec_0 = [0.98; 0; 0; 0; 1.7; 0];
tSpan = [0 8];
[t_ii, X_ii] = integrateCr3bp(xVec_0, mu_earthLuna, tSpan, tolerancePower);

figure; plot3(X_ii(:, 1), X_ii(:, 2), X_ii(:, 3))

function [t, X] = integrateCr3bp(xVec_0, mu, tSpan, tolerancePower)
    solverOptions = odeset('stats', 'off', 'RelTol', 10^tolerancePower, 'AbsTol', 10^tolerancePower);
    [t, X] = ode45(@(t, X) CR3BPUtilities.cr3bpEom(t, X, mu), tSpan, xVec_0, solverOptions);
end

function plotCr3bp(X, mu)
    figure
    hold on
    plot3(X(:, 1), X(:, 2), X(:, 3), 'DisplayName', 'Trajectory')
    plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'DisplayName', 'Earth');
    plot3(1 - mu, 0, 0, 'ro', 'MarkerSize', 5, 'DisplayName', 'Luna');
    legend;
    view(3)
    % axis equal
    hold off
end