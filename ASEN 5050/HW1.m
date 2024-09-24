mu_sun=1.32712428e11;

r=1.6599e8;
a=2.3132e8;
v_r=-11.6485;

format shortg

% 1a
nrg_mech=astroUtilities.specificNRGFromA(a, mu_sun)
v=astroUtilities.velocityAtRFromSpecificNRG(nrg_mech, mu_sun, r)

% 1b
v_theta=sqrt(v^2 - v_r^2)

% 1c
h = cross([r 0 0], [v_r, v_theta 0])
p = astroUtilities.semiLatusRectumFromH(norm(h), mu_sun)
e = astroUtilities.eccentricityFromSemiLatusRectum(p, a)

% 1d
R_1d = [1.0751e8 -1.2647e8 1.3644e5]; % km
V_1d = [1.5180e1 2.8193e1 1.0504e-2]; % km/s

H_1d = cross(R_1d, V_1d)
h_1d = norm(H_1d)
E_1d = astroUtilities.EccentricityFromRVH(R_1d, V_1d, H_1d, mu_sun)
e_1d = norm(E_1d)
nrg_mech_1d = astroUtilities.specificNRGFromRV(R_1d, V_1d, mu_sun)

% 2a
mu_jupiter=1.268e8;

v_inf = 10.7527;
theta_inf=deg2rad(139.3724);

nrg_mech_2a = v_inf^2/2
a_2a = astroUtilities.semiMajorAxisFromspecificNRG(nrg_mech_2a, mu_jupiter)
e_2a = astroUtilities.eccentricityFromNuInfinity(theta_inf)
turningAngle = astroUtilities.turningAngleFromEccentricity(e_2a)
rad2deg(turningAngle)

%2b
r_p = astroUtilities.periapsisFromae(a_2a, e_2a)
v_p = astroUtilities.velocityAtRFromSpecificNRG(nrg_mech_2a, r_p, mu_jupiter)
