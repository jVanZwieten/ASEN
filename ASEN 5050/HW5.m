addpath('..')
clear

AU = astroUtilities;

%% 1
mu_moon = AU.gravityParameter_moon; % km^3/s^2

a_1 = 7045; % km
e_1 = 0.23;
nu_1 = deg2rad(-142) % rad

%% 1a
r_1 = AU.conicEquationa(a_1, e_1, nu_1)
v_1 = AU.velocityAtRFromA(r_1, a_1, mu_moon)
h_1 = AU.angularMomentumFromae(a_1, e_1, mu_moon)
phi_1 = -AU.flightPathAngleFromhrv(h_1, r_1, v_1)
vVRth_1 = AU.velocityVRthFromVPhi(v_1, phi_1)

%% 1b
deltaVVRth = [.3; -.1; 0];
vVRth_2 = vVRth_1 + deltaVVRth

%% 1d
v_2 = norm(vVRth_2)
a_2 = AU.semiMajorAxisFromRV(r_1, v_2, mu_moon)
rVRth_1 = [r_1; 0; 0];
hVRth_2 = cross(rVRth_1, vVRth_2)
h_2 = norm(hVRth_2)
e_2 = AU.eccentricityFromAH(a_2, h_2, mu_moon)
nu_2 = AU.trueAnomalyFromaer(a_2, e_2, r_1)

%% 1e
deltaomega = nu_1 - nu_2 + 2*pi

%% 1f
m_i = 1224; % kg
I_sp = 212; % s
deltaV = norm(deltaVVRth)
m_p = AU.propellantForManeuver(deltaV*1000, m_i, I_sp)

%% 2
mu_mars = AU.gravityParameter_mars; % km^3/s^2
r_1 = 6500; % km
r_1p = 5915; % km

%% 2a
a_1 = r_1;
v_1 = sqrt(mu_mars/r_1)
e_1 = AU.eccentricityFromARp(a_1, r_1p)
h_1 = AU.angularMomentumFromae(a_1, e_1, mu_mars)
phi_1 = AU.flightPathAngleFromhrv(h_1, r_1, v_1)
vVRth_1 = AU.velocityVRthFromVPhi(v_1, phi_1)

%% 2a check
r_1a = 2*a_1 - r_1p;
e_1check = AU.eccentricityFromRaRp(r_1a, r_1p)

%% 2b
r_a2 = 7888; % km
r_p2 = 5712; % km
a_2 = AU.semiMajorAxisFromRaRp(r_a2, r_p2)
v_2 = AU.velocityAtRFromA(r_1, a_2, mu_mars)
e_2 = AU.eccentricityFromRaRp(r_a2, r_p2)
h_2 = AU.angularMomentumFromae(a_2, e_2, mu_mars)
phi_2 = AU.flightPathAngleFromhrv(h_2, r_1, v_2)
dV = AU.changeInVelocity(v_1, v_2, phi_1, phi_2)

%% 2b check
vVRth_2 = AU.velocityVRthFromVPhi(v_2, phi_2);
deltaVVRth = vVRth_2 - vVRth_1;
norm(deltaVVRth)

%% 3
mu_sun = AU.gravityParameter_sun;
a_earth = AU.semiMajorAxis_earthAu*AU.astronomicalUnitToKm
a_saturn = AU.semiMajorAxis_saturnAu*AU.astronomicalUnitToKm

%% 3a
v_earth = AU.velocityCircular(a_earth, mu_sun)
v_saturn = AU.velocityCircular(a_saturn, mu_sun)
r_ax = a_saturn;
r_px = a_earth;
a_x = AU.semiMajorAxisFromRaRp(r_ax, r_px)
v_1x = AU.velocityAtRFromA(r_px, a_x, mu_sun)
v_2x = AU.velocityAtRFromA(r_ax, a_x, mu_sun)

dV_1 = v_1x - v_earth
dV_2 = v_saturn - v_2x
dV_total = dV_1 + dV_2

TOFx = AU.timeOfFlightXfer(a_x, mu_sun)
Utilities.s2y(TOFx)

%% 3b
n_saturn = AU.meanMotionFroma(a_saturn, mu_sun)
alpha = n_saturn*TOFx
phi = pi - alpha

%% 3c
r_2 = AU.astronomicalUnitToKm*11
a_x1 = AU.semiMajorAxisFromRaRp(r_2, a_earth)
v_1x1 = AU.velocityAtRFromA(a_earth, a_x1, mu_sun)
v_2x1 = AU.velocityAtRFromA(r_2, a_x1, mu_sun)
a_x2 = AU.semiMajorAxisFromRaRp(a_saturn, r_2)
v_2x2 = AU.velocityAtRFromA(r_2, a_x2, mu_sun)
v_3x2 = AU.velocityAtRFromA(a_saturn, a_x2, mu_sun)
dV_1 = v_1x1 - v_earth
dV_2 = v_2x2 - v_2x1
dV_3 = v_saturn - v_3x2
dV_total = sum(abs([dV_1 dV_2 dV_3]))

TOF_t = AU.timeOfFlightXfer(a_x1, mu_sun) + AU.timeOfFlightXfer(a_x2, mu_sun)
Utilities.s2y(TOF_t)

ratio_a = a_saturn/a_earth

%% 4c
format long
TOF_12 = 8807.81 - 7464.33

e = 0.46;
E_1 = deg2rad(333.1741577300176);
E_2 = deg2rad(26.56519924405643);

M_1 = AU.meanAnomalyFromE(E_1, e)
M_2 = AU.meanAnomalyFromE(E_2, e)
M_12 = M_2 - M_1 +2*pi
n = deg2rad(0.02211426535691544) % rad/s
TOF_12check = M_12/n