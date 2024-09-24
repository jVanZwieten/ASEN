clear
addpath('..')

%% 1
mu_saturn = 3.794e7; % km^3/s^2
r_saturn = 60268; % km

R_1 = [-7.2; 6.7; 3.1]*10^5; % km
V_1 = [2.160; -3.360; .62]; % km/s

%% 1a
r_1 = norm(R_1)
v_1 = norm(V_1)

specNRG = astroUtilities.specificNRGFromrv(r_1, v_1, mu_saturn)

r_2 = r_saturn;
R_2_rFrame = [r_2; 0; 0];
v_2 = astroUtilities.velocityAtRFromSpecificNRG(specNRG, r_2, mu_saturn)

H = astroUtilities.AngularMomentum(R_1, V_1)
h = norm(H)
phi = -astroUtilities.flightPathAngleFromhrv(h, r_2, v_2)

V_2_rFrame = v_2*[sin(phi); cos(phi); 0]

N = astroUtilities.LineOfNodesFromH(H)
n = norm(N)
E = astroUtilities.EccentricityFromRV(R_1, V_1, mu_saturn)
e = norm(E)

i = astroUtilities.inclinationFromH(H)
Omega = astroUtilities.RAANFromN(N)
omega = -astroUtilities.argumentOfPeriapsisFromNE(N, E)
nu = -astroUtilities.trueAnomalyFromhre(h, r_2, e, mu_saturn)
theta = omega + nu

C = astroUtilities.DirectionCosineMatrix(Omega, omega, nu, i)

R_2 = C*R_2_rFrame
V_2 = C*V_2_rFrame

%% 2
clear
a = 6463.8; % km
e = .45454;
r_mars = 3397.2; % km
mu_mars = 4.28284e4; % km^3/s^2

%% 2b
P = astroUtilities.period(a, mu_mars)
r_p = astroUtilities.periapsisFromae(a, e)
alt_p = r_p - r_mars

%% 2c
i = deg2rad(74.924)
Omega = deg2rad(1.241)
omega = -2*pi + deg2rad(353.31)
nu_0 = -2*pi + deg2rad(199.38)

r_0 = astroUtilities.conicEquation(a, e, nu_0)

specNRG = astroUtilities.specificNRGFromA(a, mu_mars)
v_0 = astroUtilities.velocityAtRFromSpecificNRG(specNRG, r_0, mu_mars)

h = astroUtilities.angularMomentumFromae(a, e, mu_mars)
phi = -astroUtilities.flightPathAngleFromhrv(h, r_0, v_0)

R_0rFrame = [r_0; 0; 0];
V_0rFrame = v_0*[sin(phi); cos(phi); 0]

theta = omega + nu_0
C = astroUtilities.DirectionCosineMatrix(Omega, omega, nu_0, i)

R_0 = C*R_0rFrame
V_0 = C*V_0rFrame