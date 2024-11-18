close all; clear; clc
addpath("..")
AU = astroUtilities;

%% 1
mu_moon = CelestialParameters.gravityParameter_moon;
r_ai = 5400; % km
r_pi = 1850; % km

r_af = r_ai;
r_pf = 3000; % km

r = 3500; % km

%% 1a
a_i = AU.semiMajorAxisFromRaRp(r_ai, r_pi)
e_i = AU.eccentricityFromRaRp(r_ai, r_pi)
h_i = AU.angularMomentumFromae(a_i, e_i, mu_moon)
NRG_i = AU.specificNRGFromA(a_i, mu_moon)
v_i = AU.velocityAtRFromSpecificNRG(NRG_i, r, mu_moon)
v_iCheck = AU.velocityAtRFromA(r, a_i, mu_moon)
phi_i = -AU.flightPathAngleFromhrv(h_i, r, v_i)

a_f = AU.semiMajorAxisFromRaRp(r_af, r_pf)
e_f = AU.eccentricityFromRaRp(r_af, r_pf)
h_f = AU.angularMomentumFromae(a_f, e_f, mu_moon)
NRG_f = AU.specificNRGFromA(a_f, mu_moon)
v_f = AU.velocityAtRFromSpecificNRG(NRG_f, r, mu_moon)
phi_f = AU.flightPathAngleFromhrv(h_f, r, v_f)

deltaPhi1 = phi_f - phi_i
deltaPhi2 = -phi_f - phi_i

deltaV1 = AU.changeInVelocityDPhi(v_i, v_f, deltaPhi1)
deltaV2 = AU.changeInVelocityDPhi(v_i, v_f, deltaPhi2)

%% 1b
theta_i = -AU.trueAnomalyFromaer(a_i, e_i, r)
theta_f = -AU.trueAnomalyFromaer(a_f, e_f, r)
deltaomega = theta_i - theta_f

%% 2
mu_jupiter = CelestialParameters.gravityParameter_jupiter;
m_i = 2500; % kg
m_pi = 500; % kg
I_sp = 300; % kg
mu_X = 4000; % km^3/s^2

%% 2a
vVecRt_in = [2.2789; 5.8841]; % km/s
a_in = 1.8e6; % km

v_in = norm(vVecRt_in)
NRG_in = AU.specificNRGFromA(a_in, mu_jupiter)
r_in = AU.radialDistanceAtVFromSpecificNRG(v_in, NRG_in, mu_jupiter)

%% 2b
v_X = AU.velocityCircular(r_in, mu_jupiter)
vVecRt_X = [0; v_X];
vVecRt_infIn = vVecRt_in - vVecRt_X

%% 2c
v_out = 7.6759; % km/s
phi_out = deg2rad(20.91);
vVecRt_out = v_out*[sin(phi_out); cos(phi_out)]

%% 2d
vVecRt_infOut = vVecRt_out - vVecRt_X
Utilities.assertIsWithin(norm(vVecRt_infOut), norm(vVecRt_infIn), -3)

%% 2h
v_inf = norm(vVecRt_infIn)
vHat_infIn = vVecRt_infIn/v_inf
vHat_infOut = vVecRt_infOut/v_inf
del = AU.turningAngleFromVVecInfInOut(vVecRt_infIn, vVecRt_infOut)
Utilities.assertIsWithin(del, AU.turningAngleFromVVecInfInOut(vHat_infIn, vHat_infOut), -5)
e_h = AU.eccentricityFromTurningAngle(del)

a_h = AU.semiMajorAxisHyperbolaFromVInfinity(v_inf, mu_X)
r_ph = AU.periapsisFromae(a_h, e_h)

%% 2i
dVVecRt = vVecRt_out - vVecRt_in
dV = norm(dVVecRt)

%% 2j
dV_max = AU.idealRocketEquation(I_sp, m_i, m_i - m_pi)

%% 3
mu_sun = CelestialParameters.gravityParameter_sun;
r_jupiter = CelestialParameters.semiMajorAxis_jupiter

vVecRt_in = [3.2476; 13.1175]; % km/s
vVecRt_out = [3.3426; 13.7356]; % km/s

v_jupiter = AU.velocityCircular(r_jupiter, mu_sun)
vVecRt_jupiter = [0; v_jupiter];
vVecRt_infIn = vVecRt_in - vVecRt_jupiter
vVecRt_infOut = vVecRt_out - vVecRt_jupiter

v_infIn = norm(vVecRt_infIn)
v_infOut = norm(vVecRt_infOut)
Utilities.assertIsWithin(v_infIn, v_infOut, -5)