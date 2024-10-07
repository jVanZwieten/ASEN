clear
addpath("..")

mu_mars = 4.305e4; % km^3/s^2
r_mars = 3397.2; % km

%% 1
h = 2.4799e4; % km^2/s
r_1 = 15000; % km
v_1 = 1.8111; % km/s
r_2 = 19000; % km

%% 1a
specNRG = astroUtilities.specificNRGFromrv(r_1, v_1, mu_mars)
a = astroUtilities.semiMajorAxisFromspecificNRG(specNRG, mu_mars)

%% 1b
e = astroUtilities.eccentricityFromspecNRGh(specNRG, h, mu_mars)

%% 1c
nu_1 = -astroUtilities.trueAnomalyFromhre(h, r_1, e, mu_mars)
E_1 = astroUtilities.eccentricAnomalyFromnu(nu_1, e)

% as a check:
E_1check = -astroUtilities.eccentricAnomalyFromaer(a, e, r_1)
n_1check = astroUtilities.trueFromEccentricAnomaly(E_1check, e)
n_1check0 = -astroUtilities.trueAnomalyFromaer(a, e, r_1)

%% 1d
nu_2 = -astroUtilities.trueAnomalyFromaer(a, e, r_2)
nu_2check = -astroUtilities.trueAnomalyFromhre(h, r_2, e, mu_mars)
E_2 = astroUtilities.eccentricAnomalyFromnu(nu_2, e)

%% 1e
t_12 = astroUtilities.timeBetweenEccentricAnomaliesa(E_1, E_2, a, e, 1, mu_mars)
T = astroUtilities.period(a, mu_mars)
Utilities.s2days(t_12)

rad2deg(nu_1)
rad2deg(E_1)
rad2deg(nu_2)
rad2deg(E_2)

%% 2
rV_3 = [-7.6650e3; 6.5468e3; -4.5740e2]; % km
vV_3 = [1.6334; 0.1226; -1.9455]; % km/s

%% 2a
hV = astroUtilities.AngularMomentum(rV_3, vV_3)
h = norm(hV)
nV = astroUtilities.LineOfNodesFromH(hV)
n = norm(nV)
r_3 = norm(rV_3)
v_3 = norm(vV_3)
eV = astroUtilities.EccentricityFromRV(rV_3, vV_3, mu_mars)
e = norm(eV)
rdotv3 = dot(rV_3, vV_3)

a = astroUtilities.semiMajorAxisFromhe(h, e, mu_mars)
i = astroUtilities.inclinationFromH(hV)
Omega = -astroUtilities.RAANFromN(nV)
omega = -astroUtilities.argumentOfPeriapsisFromNE(nV, eV)
nu_3 = -astroUtilities.trueAnomalyFromRE(rV_3, eV)

r_3check = astroUtilities.conicEquationa(a, e, nu_3);
assert(within(r_3, r_3check, -3))
rVRnh_3 = [r_3check; 0; 0];
specNRGcheck = astroUtilities.specificNRGFromA(a, mu_mars);
v_3check = astroUtilities.velocityAtRFromSpecificNRG(specNRGcheck, r_3check, mu_mars);
assert(within(v_3, v_3check, -3))
phi_3 = -astroUtilities.flightPathAngleFromhrv(h, r_3check, v_3check);
vVRnh_3 = [sin(phi_3); cos(phi_3); 0]*v_3check;
C = astroUtilities.DirectionCosineMatrix(Omega, omega, nu_3, i);

rV_3check = C*rVRnh_3
vV_3check = C*vVRnh_3

%% 2b
t_34 = Utilities.h2s(2) % s
meanMotion = astroUtilities.meanMotionFroma(a, mu_mars)
E_3 = astroUtilities.eccentricAnomalyFromnu(nu_3, e)
M_3 = astroUtilities.meanAnomalyFromE(E_3, e)

M_4 = M_3 + meanMotion*t_34
E_4 = astroUtilities.KeplerNewtonSolver(M_4, a, e, mu_mars)
nu_4 = astroUtilities.trueFromEccentricAnomaly(E_4, e)

p = astroUtilities.semiLatusRectumFromae(a, e)
dnu_34 = nu_4 - nu_3
r_4 = astroUtilities.conicEquationp(p, e, nu_4)

f = astroUtilities.fFunction(r_4, p, dnu_34)
fDot = astroUtilities.fDotFunction(mu_mars, p, dnu_34, r_4, r_3)
g = astroUtilities.gFunction(r_4, r_3, mu_mars, p, dnu_34)
gDot = astroUtilities.gDotFunction(r_3, p, dnu_34)

rV_4 = f*rV_3 + g*vV_3
vV_4 = fDot*rV_3 + gDot*vV_3

assert(within(r_4, norm(rV_4), -3))
rVRnh_4 = [r_4; 0; 0];
specNRGcheck = astroUtilities.specificNRGFromA(a, mu_mars);
v_4 = astroUtilities.velocityAtRFromSpecificNRG(specNRGcheck, r_4, mu_mars);
assert(within(v_4, norm(vV_4), -3))
phi_4 = astroUtilities.flightPathAngleFromhrv(h, r_4, v_4);
vVRnh_4 = [sin(phi_4); cos(phi_4); 0]*v_4;
C = astroUtilities.DirectionCosineMatrix(Omega, omega, nu_4, i);

rV_4check = C*rVRnh_4
vV_4check = C*vVRnh_4

%% 3
rHat_dn = [-.64279; -.76604; 0];
rHat_a = [-.02970; -.97508; -.21985];

%% 3a
nHat = -rHat_dn
Omega = acos(nHat(1))

%% 3b
eHat = -rHat_a
omega = acos(dot(nHat, eHat))

function result = within(a, b, d)
    result = abs(a - b) < 10^d;
end