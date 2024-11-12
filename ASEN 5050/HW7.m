clear; clc;
addpath("..")
AU = astroUtilities;

%% 1
mu_saturn = CelestialParameters.gravityParameter_saturn;
r_p = 6e5; % km
r_a = 1.8e6; % km
a_titan = 1.22183e6; % km
m_titan = 1.3455e23; % kg
R_titan = 2575; % km

%% 1a
r_1 = a_titan;
a = AU.semiMajorAxisFromRaRp(r_a, r_p)
e = AU.eccentricityFromRaRp(r_a, r_p)
nu_1 = -AU.trueAnomalyFromaer(a, e, r_1)
v_1 = AU.velocityAtRFromA(r_1, a, mu_saturn)
phi_1 = AU.flightPathAngleFromENu(e, nu_1)
vVecRth_1 = v_1*[sin(phi_1); cos(phi_1); 0]

v_titan = AU.velocityCircular(a_titan, mu_saturn)
vVecRth_titan = [0; v_titan; 0]
vVecRth_infIn = vVecRth_1 - vVecRth_titan

%% 1b
mu_titan = CelestialParameters.gravitationalParameter(m_titan)
r_pH = 3000; % km
v_infIn = norm(vVecRth_infIn)
a_h = AU.semiMajorAxisHyperbolaFromVInfinity(v_infIn, mu_titan)
e_h = AU.eccentricityFromARp(a_h, r_pH)
del = AU.turningAngleFromEccentricity(e_h)

%% 1e
vVecRt_in = vVecRth_1(1:2, 1)
vVecRt_infIn = vVecRth_infIn(1:2, 1)
vVecRt_infOut = [cos(del) -sin(del); sin(del) cos(del)]*vVecRt_infIn
vVecRt_out = vVecRt_in - vVecRt_infIn + vVecRt_infOut
vVecRt_outCheck = AU.velocityAfterFlyby(vVecRt_in, vVecRth_titan(1:2, 1), del)

%% 1f
rVecRth_out = [r_1; 0; 0];
vVecRth_out = [vVecRt_out; 0];
eVecRth_out = AU.EccentricityFromRV(rVecRth_out, vVecRth_out, mu_saturn)
e_out = norm(eVecRth_out)
v_out = norm(vVecRth_out)
a_out = AU.semiMajorAxisFromRV(r_1, v_out, mu_saturn)
nu_out = -AU.trueAnomalyFromaer(a_out, e_out, r_1)

%% 1g
deltaVVecRth = vVecRth_out - vVecRth_1

%% 2e
r_1 = [2.986942867720575E+08; 2.010772111651812E+08; -4.192849591564223E+06];
v_1 = [5.289894007442628; 2.209289843371483E+01; -2.627701845286099E-01];
r_f = [2.98703e+08; 2.01087e+08; -4.19285e+06];
v_f = [5.28979; 22.0977; -0.26277];

deltaR = r_f - r_1
deltaV = v_f - v_1