close all; clear; clc; addpath("../")
AU = astroUtilities;

mu_earth = CelestialParameters.gravityParameter_earth;
a = 10000; % km

%% 2
rVec_0 = [2; 2];
vVec_0minus = [-.03; .01; .05];
rVec_1 = [-2; 2];
vVec_1plus = [0; 0; 0];

%% 2a
n = AU.meanMotion(a, mu_earth)

Phi_rr1 = [7 0; -6*pi 1]
Phi_rv1 = [0 4/n; -4/n, -3*pi/n]

nt_1 = pi;
Phi_rrCheck = [4-3*cos(nt_1) 0; 6*(sin(nt_1)-nt_1) 1];
Phi_rvCheck = [1/n*sin(nt_1) 2/n*(1-cos(nt_1)); 2/n*(cos(nt_1)-1) 1/n*(4*sin(nt_1)-3*nt_1)];

Utilities.assertIsWithin(Phi_rr1, Phi_rrCheck, 1e-6);
Utilities.assertIsWithin(Phi_rv1, Phi_rvCheck, 1e-6);

Phi_rv1inv = inv(Phi_rv1)
vVec_0plus = Phi_rv1inv*(rVec_1 - Phi_rr1*rVec_0)
vVec_0PlusCheck = Phi_rv1\(rVec_1 - Phi_rrCheck*rVec_0);
Utilities.assertIsWithin(vVec_0plus, vVec_0PlusCheck, 1e-6);

vVec_0plus = [vVec_0plus; 0]
dVVec_0 = vVec_0plus - vVec_0minus
dV_0 = norm(dVVec_0)

%% 2b
Phi_vr1 = [3*n*sin(nt_1) 0 0; 6*n*(cos(nt_1)-1) 0 0; 0 0 -n*sin(nt_1)]
Phi_vv1 = [cos(nt_1) 2*sin(nt_1) 0; -2*sin(nt_1) 4*cos(nt_1)-3 0; 0 0 cos(nt_1)]
rVec_0 = [rVec_0; 0]; % m

vVec_1minus = Phi_vr1*rVec_0 + Phi_vv1*vVec_0plus
dVVec_1 = vVec_1plus - vVec_1minus
dV_1 = norm(dVVec_1)

%% 2c
rVec_0 = rVec_0(1:2);
vVec_0plus = vVec_0plus(1:2);

nt = [0:.1:pi, pi];
X = zeros(2, length(nt));
for i = 1:length(nt)
    nt_i = nt(i);
    Phi_rr = [4-3*cos(nt_i) 0; 6*(sin(nt_i)-nt_i) 1];
    Phi_rv = [1/n*sin(nt_i) 2/n*(1-cos(nt_i)); 2/n*(cos(nt_i)-1) 1/n*(4*sin(nt_i)-3*nt_i)];
    X(:, i) = Phi_rr*rVec_0 + Phi_rv*vVec_0plus;
end

figure
plot(X(2, :), X(1, :)); xlabel("y (km)"); ylabel("x (km)"); title("Relative Motion of CubeSat in LVLH Frame"); axis equal
hold on
plot(X(2, 1), X(1, 1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(X(2, end), X(1, end), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
legend('Trajectory', 'Initial Position', 'Final Position');
ylim([-2.5 2.5]);


%% 3
P_rotMars = CelestialParameters.rotationPeriod_mars
mu_mars = CelestialParameters.gravityParameter_mars;

omega_rotMars = CelestialParameters.rotationRate_mars
dLambda = 5*pi/6;
T = dLambda/omega_rotMars
a = AU.semiMajorAxisFromPeriod(T, mu_mars)