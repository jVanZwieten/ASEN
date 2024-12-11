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
n = AU.meanMotion(mu_earth, a)

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

vVec_1minus = Phi_vr1*rVec_0 + Phi_vv1*vVec_0minus
dVVec_1 = vVec_1plus - vVec_1minus
dV_1 = norm(dVVec_1)

%% 3
P_rotMars = CelestialParameters.rotationPeriod_mars
mu_mars = CelestialParameters.gravityParameter_mars;

omega_rotMars = CelestialParameters.rotationRate_mars
dLambda = 5*pi/6;
T = dLambda/omega_rotMars
a = AU.semiMajorAxisFromPeriod(mu_mars, T)