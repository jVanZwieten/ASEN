clear; clc;
addpath("..")

AU = astroUtilities;

mu_sun = CelestialParameters.gravityParameter_sun;
mu_mars = CelestialParameters.gravityParameter_mars;

r_c = CelestialParameters.au2km(1.52);
v_c = AU.velocityCircular(r_c, mu_sun);
v_cExpected = 24.1586; % km/s
Utilities.assertIsWithin(v_c, v_cExpected, -4)

a_x = CelestialParameters.au2km(5.75); % km
e_x = 0.8104;
nu = AU.trueAnomalyFromaer(a_x, e_x, r_c);
nu_expected = deg2rad(68.4);
Utilities.assertIsWithin(nu, nu_expected, -3)

h_x = AU.angularMomentumFromae(a_x, e_x, mu_sun);
v_r = AU.velocityRadial(h_x, e_x, nu, mu_sun);
v_rExpected = 15.97; %km/s
Utilities.assertIsWithin(v_r, v_rExpected, -2)

v_theta = AU.velocityTangential(h_x, e_x, nu, mu_sun);
v_thetaExpected = 27.53; % km/s
Utilities.assertIsWithin(v_theta, v_thetaExpected, -2)

vVec_in = [v_r; v_theta];
vVec_mars = [0; v_c];
vVec_infIn = vVec_in - vVec_mars;

v_inf = norm(vVec_infIn);
a_h = AU.semiMajorAxisHyperbolaFromVInfinity(v_inf, mu_mars);
a_hExpected = -161.54; % km
Utilities.assertIsWithin(a_h, a_hExpected, -2)

r_p = 200 + CelestialParameters.radius_Mars;
e_h = astroUtilities.eccentricityFromARp(a_h, r_p);
e_hExpected = 23.27;
Utilities.assertIsWithin(e_h, e_hExpected, -2)

del = AU.turningAngleFromEccentricity(e_h);
delExpected = deg2rad(4.926);
Utilities.assertIsWithin(del, delExpected, -3)

vVec_out = astroUtilities.velocityAfterFlyby(vVec_in, vVec_mars, -del);
v_out = norm(vVec_out);
v_outExpected = 30.76; % km/s
Utilities.assertIsWithin(v_out, v_outExpected, -2);