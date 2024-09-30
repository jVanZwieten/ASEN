clear
addpath('..')

%% 1
mu_saturn = 3.794e7; % km^3/s^2
r_saturn = 60268; % km

R_1 = [-7.2; 6.7; 3.1]*10^5; % km
V_1 = [2.160; -3.360; .62]; % km/s

%% 1a
% t_1 = initial epoch
% t_i = impact epoch
r_1 = norm(R_1)
v_1 = norm(V_1)

r_i = r_saturn;

H = astroUtilities.AngularMomentum(R_1, V_1);
h = norm(H);
e = norm(astroUtilities.EccentricityFromRVH(R_1, V_1, H, mu_saturn))

p = astroUtilities.semiLatusRectumFromh(h, mu_saturn)
dot(R_1, V_1)
nu_1 = -astroUtilities.trueAnomalyFrompre(p, r_1, e)
nu_i = -astroUtilities.trueAnomalyFrompre(p, r_i, e)
dnu = nu_i - nu_1

f = astroUtilities.fFunction(r_i, p, dnu)
fDot = astroUtilities.fDotFunction(mu_saturn, p, dnu, r_i, r_1)
g = astroUtilities.gFunction(r_i, r_1, mu_saturn, p, dnu)
gDot = astroUtilities.gDotFunction(r_1, p, dnu)

R_i = astroUtilities.positionfg(R_1, V_1, f, g)
V_i = astroUtilities.velocityfDotgDot(R_1, V_1, fDot, gDot)

%% 1b
a = astroUtilities.semiMajorAxisFrompe(p, e)
n = astroUtilities.meanMotionFroma(a, mu_saturn)

E_1 = astroUtilities.eccentricAnomalyFromnu(nu_1, e)
E_i = astroUtilities.eccentricAnomalyFromnu(nu_i, e)

t_1i = astroUtilities.timeBetweenEccentricAnomalies(E_1, E_i, e, n)
t_1i/60/60


%% 2
mu_jupiter = 1.268e8; % km^3/s^2
r_jupiter = 71492; % km
R_1 = [5.352950e6; 7.053778e5; -4.0597e5]; % km
V_1 = [-4.164248; 1.963690; 3.191257e-1]; % km/s

%% 2a
r_1 = norm(R_1)
v_1 = norm(V_1)

E = astroUtilities.EccentricityFromRV(R_1, V_1, mu_jupiter)
e = norm(E)

dot(R_1, V_1)
nu_1 = -astroUtilities.trueAnomalyFromRE(R_1, E)
E_1 = astroUtilities.eccentricAnomalyFromnu(nu_1, e)

%% 2b
H = astroUtilities.AngularMomentum(R_1, V_1)
N = astroUtilities.LineOfNodesFromH(H)
n = norm(N)

omega = astroUtilities.argumentOfPeriapsisFromNE(N, E)
nu_2 = pi - omega

E_2 = astroUtilities.eccentricAnomalyFromnu(nu_2, e)

%% 2c
specNRG = astroUtilities.specificNRGFromrv(r_1, v_1, mu_jupiter)
a = astroUtilities.semiMajorAxisFromspecificNRG(specNRG, mu_jupiter)
meanMotion = astroUtilities.meanMotionFroma(a, mu_jupiter)

t_12 = astroUtilities.timeBetweenEccentricAnomalies(E_1, E_2, e, meanMotion)
t_12_days = Utilities.s2days(t_12)

%% 2d
h = norm(H)
p = astroUtilities.semiLatusRectumFromh(h, mu_jupiter)
deltaNu = nu_2 - nu_1
r_2 = astroUtilities.conicEquationp(p, e, nu_2)

f = astroUtilities.fFunction(r_2, p, deltaNu)
fDot = astroUtilities.fDotFunction(mu_jupiter, p, deltaNu, r_2, r_1)
g = astroUtilities.gFunction(r_2, r_1, mu_jupiter, p, deltaNu)
gDot = astroUtilities.gDotFunction(r_1, p, deltaNu)

R_2 = f*R_1 + g*V_1
V_2 = fDot*R_1 + gDot*V_1