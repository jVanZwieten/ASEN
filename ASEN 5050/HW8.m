clear; clc; close all; addpath("..");
AU = astroUtilities;

%% 1
mu_sun = CelestialParameters.gravityParameter_sun;

r_p = 7500; % km
r_a = 8500; % km
i = deg2rad(105); % rad
T = 110*60 % s
R = 6500; % km
r_c = 2.25; % AU

a = AU.semiMajorAxisFromRaRp(r_a, r_p)
mu = 4*pi^2*a^3/T^2
G = CelestialParameters.universalGravitationalConstant
M = mu/G
r_c = CelestialParameters.au2km(r_c) % km
T = AU.period(r_c, mu_sun)
OmegaDot = 2*pi/T
e = AU.eccentricityFromRaRp(r_a, r_p)
J_2 = AU.J2FromOmegaDot(OmegaDot, R, mu, a, e, i)

%% 2
RVec_0 = [2489.63813; -3916.07418; -5679.05521]; % km
VVec_0 = [9.13452; -1.91212; 2.57306]; % km/s
XVec_0 = [RVec_0; VVec_0];

%% 2a
mu_earth = CelestialParameters.gravityParameter_earth;
R_0 = norm(RVec_0)
V_0 = norm(VVec_0)
nrg_0 = AU.specificNRGFromRV(R_0, V_0, mu_earth)
hVec_0 = cross(RVec_0, VVec_0)
h_0 = norm(hVec_0)

%% 2b
p = AU.semiLatusRectumFromh(h_0, mu_earth)
e = AU.eccentricityFromspecNRGh(nrg_0, h_0, mu_earth)
dot(RVec_0, VVec_0)
nu_0 = AU.trueAnomalyFrompre(p, R_0, e)
nu_1 = pi;
deltaNu = nu_1 - nu_0
R_1 = AU.conicEquationp(p, e, nu_1)
f = AU.fFunction(R_1, p, deltaNu)
fDot = AU.fDotFunction(mu_earth, p, deltaNu, R_1, R_0)
g = AU.gFunction(R_1, R_0, mu_earth, p, deltaNu)
gDot = AU.gDotFunction(R_0, p, deltaNu)

C = [f*eye(3) g*eye(3); fDot*eye(3) gDot*eye(3)];
XVec_ref = C*XVec_0

RVec_ref = f*RVec_0 + g*VVec_0
VVec_ref = fDot*RVec_0 + gDot*VVec_0

E_0 = AU.eccentricFromTrueAnomaly(nu_0, e)
E_1 = AU.eccentricFromTrueAnomaly(nu_1, e)

M_0 = AU.meanAnomalyFromE(E_0, e)
M_1 = AU.meanAnomalyFromE(E_1, e)
deltaM = M_1 - M_0
a = AU.semiMajorAxisFrompe(p, e)
n = AU.meanMotionFroma(a, mu_earth)
deltaT = deltaM/n

deltaTCheck = AU.timeBetweenEccentricAnomalies(E_0, E_1, e, n, 0)
Utilities.assertIsWithin(deltaT, deltaTCheck, 1e-10)

%% 2c
tspan = [0 deltaT];
solverOptions = odeset('stats', 'off', 'RelTol', 1e-8, 'AbsTol', 1e-8);
[t, X] = ode45(@(t, X) twoBodyODE(t, X, mu_earth), tspan, XVec_0, solverOptions);

RM_numerical = X(:, 1:3);
VM_numerical = X(:, 4:6);

figure;
plot3(RM_numerical(:,1), RM_numerical(:,2), RM_numerical(:,3));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('Spacecraft Trajectory');
grid on;
hold on;
plot3(XVec_0(1), XVec_0(2), XVec_0(3), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
plot3(0, 0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
midT = round(length(t)/2);
quiver3(RM_numerical([1 midT], 1), RM_numerical([1 midT], 2), RM_numerical([1 midT], 3), VM_numerical([1 midT], 1), VM_numerical([1 midT], 2), VM_numerical([1 midT], 3), 0.5, 'r');
legend('Trajectory', 'Initial Position', 'Earth');
axis equal;

NrgHM = zeros(length(t), 2);
for i = 1:length(t)
    NrgHM(i, 1) = AU.specificNRGFromRV(norm(RM_numerical(i,:)), norm(VM_numerical(i,:)), mu_earth);
    NrgHM(i, 2) = norm(AU.angularMomentum(RM_numerical(i,:), VM_numerical(i,:)));
end

figure;
tiledlayout(2,1);
nexttile;
plot(t, NrgHM(:,1));
ylabel('Specific Energy (km^2/s^2)');
nexttile
plot(t, NrgHM(:,2));
ylabel('Angular Momentum (km^2/s)');
xlabel('Time (s)');


%% 2d
numericalData = table('Size', [0 6], 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'tolerance', 'deltaR', 'deltaV', 'deltaNrg', 'deltaH', 'compTime'});
for tol = -4:-2:-12
    solverOptions = odeset('stats', 'off', 'RelTol', 10^tol, 'AbsTol', 10^tol);
    tic
    [t, X] = ode45(@(t, X) twoBodyODE(t, X, mu_earth), tspan, XVec_0, solverOptions);
    compTime = toc;

    RVec_1 = X(end, 1:3)';
    deltaR = norm(RVec_ref - RVec_1);
    VVec_1 = X(end, 4:6)';
    deltaV = norm(VVec_ref - VVec_1);

    RVec_0 = X(1, 1:3)';
    VVec_0 = X(1, 4:6)';
    deltaNrg = AU.specificNRGFromRV(RVec_1, VVec_1, mu_earth) - AU.specificNRGFromRV(RVec_0, VVec_0, mu_earth);
    deltaH = norm(AU.angularMomentum(RVec_1, VVec_1)) - norm(AU.angularMomentum(RVec_0, VVec_0));
    numericalData = [numericalData; {tol, deltaR, deltaV, deltaNrg, deltaH, compTime}];
end

%% 2f
R_earth = CelestialParameters.radius_earth/1000; % km
J_2 = CelestialParameters.j2_earth;

i = deg2rad(63.4);
omega = deg2rad(270);
OmegaDot = AU.RAANDotFromJ2(J_2, R_earth, mu_earth, a, e, i)
omegaDot = AU.argumentOfPeriapsisDotFromJ2(J_2, R_earth, mu_earth, a, e, i)
OmegaDot_sunSynch = 2*pi/(365.25*24*3600)

function dXdt = twoBodyODE(t, X, mu)
    R = X(1:3);
    V = X(4:6);
    r = norm(R);
    a = -mu / r^3 * R;
    dXdt = [V; a];
end