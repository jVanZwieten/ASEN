clear; clc; close all; addpath('..')
HU = HypersonicsUtilities;

gamma = HU.ratioSpecificHeats_air;
R = HU.gasConstant_air;
c_p = HU.specificHeatAtConstantPressure_air;

%% HW6 problem 3
h = 35; % km
M_inf = 10;
T_w = 1000; % K
R_N = .01; % m

%% a
rho_sl = 1.225; % kg/m^3
rhoRatio_35 = 6.9089e-3; % kg/m^3
rho_inf = rho_sl*rhoRatio_35
T_inf = 236.513; % K

a_inf = HU.speedOfSoundAir(T_inf)
u_inf = a_inf*M_inf
p_inf = HU.idealGasPressure(rho_inf, T_inf)
mu_inf = HU.sutherlandsLawViscosity(T_inf)

%% b
atm = HU.airPressure_seaLevel;

p_2 = p_inf + rho_inf*u_inf^2
p_2atm = p_2/atm
h_1 = c_p*T_inf
h_2 = h_1 + u_inf^2/2

T_2 = 3275; % K
z_2 = 1.055;

rho_2 = HU.densityRealGasEffect(p_2, z_2, R, T_2)
u_2 = rho_inf*u_inf/rho_2

%% c
T_e = 3300; % K
z_e = 1.055;
SoverR = 35;
p_e = 8.1060e4; % Pa
rho_e = 8.1126e-2; % kg/m^3
mu_e = 8.1044e-5; % kg/m/s

%% d
p_w = p_e;
rho_w = HU.idealGasDensity(p_w, T_w)
mu_w = HU.sutherlandsLawViscosity(T_w)

%% e
cylinderCorrection = 0.747;
q_w = HU.detraEtAllHeatFlux(rho_inf, u_inf, R_N)*cylinderCorrection

%% f
F = 1;
Pr_w = 0.75;

k = 0.57;
dudx_e = HU.velocityGradient(R_N, p_e, p_inf, rho_e)
h_0e = 4.9891e6; % J/kg
h_w = c_p*T_w

q_w = HU.fayRiddellHeatFlux(k, Pr_w, rho_e, mu_e, rho_w, mu_w, dudx_e, h_0e, h_w, F)

%% g
rho_d = 150000;
theta_d = 59500
alphaOverAlpha1 = HU.massFractionFunction(rho_e, rho_d, theta_d, T_e)
alpha = roots([1, alphaOverAlpha1, -alphaOverAlpha1])
C_allO = .22;
C_O = alpha(2)*C_allO
C_O2 = C_allO - C_O

%% midterm
%% 1
R_N = 4.8; % m
R_B = 2.25; % m
A = pi*R_B^2
m = 7500; % kg
c_D = 1.6;
L = 0;
T = 0;
k = 2e-4/sqrt(R_N)
V_0 = 11100; % m/s
h_0 = 100; % km
gamma_0 = deg2rad(-7); % rad

%% a
beta = HU.beta(m, gamma_0, c_D, A)
q_max = HU.maxHeatingBallistic(k, beta, V_0)
g_max = HU.maxG(V_0, gamma_0)

%% 2
u_inf = 5000; % m/s
h = 45; % km
T_inf = 264.164; % K
p_inf = 149.1; % Pa
rho_inf = 1.9662e-3; % kg/m^3

%% 2a
p_2 = p_inf + rho_inf*u_inf^2
p_2atm = p_2/atm
h_2 = c_p*T_inf + u_inf^2/2

T_2 = 5500; % K
z_2 = 1.27;

rho_2 = HU.densityRealGasEffect(p_2, z_2, R, T_2)
u_2 = rho_inf*u_inf/rho_2

%% 2b
mu_2 = HU.sutherlandsLawViscosity(T_2)
mu_2 = 1.1*mu_2
k_2 = HU.sutherlandsLawThermalConductivity(T_2)

%% 3c
k = 10;
R = 6.5;
r = 1.5;
m = k*4*pi^2*R*r

C_pMax = 1.839;
CDA = 8*pi/3*R*r*C_pMax
CDA_apollo = 1.6*A

%% 1
T_w = 1000; % K
Pr = 0.75;
M_inf = 10;
R_N = 0.01; % m
L = 4; % m

%% a
T_inf = 250.35; % K
rhoRatio = 3.2618e-3;
rho_inf = HU.airDensity_seaLevel*rhoRatio

a_inf = HU.speedOfSoundAir(T_inf)
u_inf = a_inf*M_inf
p_inf = HU.idealGasPressure(rho_inf, T_inf)
mu_inf = HU.sutherlandsLawViscosity(T_inf)
Re_inf = HU.reynoldsNumber(rho_inf, u_inf, L, mu_inf)

%% b
T_0e = T_inf*HU.stagnationTemperatureRatioAir(M_inf)

%% 1
T_e = 500; % K
M_e = 15;
rho_e = .5; % kg/m^3
T_w = 1000; % K
Pr_w = 0.75;
L = 5; % m

%% a
p_e = HU.idealGasPressure(rho_e, T_e)
a_e = HU.speedOfSoundAir(T_e)
u_e = a_e*M_e
mu_e = HU.sutherlandsLawViscosity(T_e)

%% b
r = sqrt(Pr)
T_0e = T_e*HU.stagnationTemperatureRatioAir(M_e)
TStar = HU.EckertReferenceTemperature(T_e, T_0e, T_w, r)
rhoStar = HU.idealGasDensity(p_e, TStar)
muStar = HU.sutherlandsLawViscosity(TStar)
ReStar = rhoStar*u_e/muStar
C_fx = sqrt(3)*0.664/sqrt(ReStar)

C_h = HU.reynoldsAnalogy(C_fx, Pr)
h_aw = HU.hawFromR(r, c_p*T_0e, c_p*T_e)
qDot_w = HU.heatFluxFromHeatingCoefficient(C_h, rhoStar, u_e, h_aw, c_p*T_w)

%% c
e = 0.95;
qDot_rad = HU.radiationHeatFlux(T_w, e)
x_eq = (qDot_w/qDot_rad)^2