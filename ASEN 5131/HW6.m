close all; clear; clc; addpath('..'); format shortG
HU = HypersonicsUtilities;

%% 3
alt = 35000; % m
M_1 = 10;
T_w = 1000; % K
R_N = .01; % m

%% 3a
rhoRatio = 6.9089e-3;
T_1 = 236.513; % K
rho_1 = HU.airDensity_seaLevel*rhoRatio

R = HU.gasConstant_air;
gamma = 1.4;
a = HU.speedOfSound(gamma, R, T_1)

u_1 = M_1*a
p_1 = rho_1*R*T_1
mu_1 = HU.sutherlandsLawViscosity(T_1)

%% 3b
p_2 = p_1 + rho_1*u_1^2
atm = HU.airPressure_seaLevel;
p_2atm = p_2/atm
c_p = HU.specificHeatAtConstantPressure(gamma, R)
h_2 = c_p*T_1 + .5*u_1^2
T_2 = 3250; % K
z_2 = 1.055;
rho_2 = HU.densityRealGasEffect(p_2, z_2, R, T_2)
u_2 = rho_1*u_1/rho_2

p_2 = p_1 + rho_1*u_1^2 - rho_2*u_2^2
p_2atm = p_2/atm
h_2 = c_p*T_1 + .5*u_1^2 - .5*u_2^2
T_2 = 3225; % K
z_2 = 1.055;
rho_2 = HU.densityRealGasEffect(p_2, z_2, R, T_2)
u_2 = rho_1*u_1/rho_2
SRratio = 35;
S = SRratio*R


%% 3c
h_0e = HU.enthalpyTotal(h_2, u_2)
p_0eAtm = 0.8; % atm
p_0e = p_0eAtm*atm
T_0e = 3275; % K
z_0e = 1.005;

rho_0e = HU.densityRealGasEffect(p_0e, z_0e, R, T_0e)
mu_0e = HU.sutherlandsLawViscosity(T_0e)

%% 3d
p_w = p_0e;
rho_w = p_w/(R*T_w)
mu_w = HU.sutherlandsLawViscosity(T_w)

%% 3e
qDot_w = HU.detraEtAllHeatFlux(rho_1, u_1, R_N)
qDot_w = .747*qDot_w

%% 3f
F = 1;
Pr_w = 0.75;
k = 0.57;
dudx_e = HU.velocityGradient(R_N, p_0e, p_1, rho_0e)
h_w = c_p*T_w

qDot_w = HU.fayRiddellHeatFlux(k, Pr_w, rho_0e, mu_0e, rho_w, mu_w, dudx_e, h_0e, h_w, F)

%% 3g
rho_dO = 1.5e5; % kg/m^3
theta_d = 5.95e4; % K
alphaOverAlpha1_O = HU.massFractionFunction(rho_0e, rho_dO, theta_d, T_0e)
alpha_O = roots([1, alphaOverAlpha1_O, -alphaOverAlpha1_O])
alpha_O = alpha_O(2)
C_O = .2*alpha_O

rho_dN = 1.3e5; % kg/m^3
theta_d = 1.13e5; % K
alphaOverAlpha1_N = HU.massFractionFunction(rho_0e, rho_dN, theta_d, T_0e)
alpha_N = roots([1, alphaOverAlpha1_N, -alphaOverAlpha1_N])
alpha_N = alpha_N(2)

h_O = 15.46e6
h_D = C_O*h_O
h_dOverH_0e = h_D/h_0e

%% 3f
Le = 1.4;
F_1 = HU.fayRiddelCase1(Le, h_dOverH_0e)
qDot_w1 = HU.fayRiddellHeatFlux(k, Pr_w, rho_0e, mu_0e, rho_w, mu_w, dudx_e, h_0e, h_w, F_1)

F_2 = HU.fayRiddelCase2(Le, h_dOverH_0e)
qDot_w2 = HU.fayRiddellHeatFlux(k, Pr_w, rho_0e, mu_0e, rho_w, mu_w, dudx_e, h_0e, h_w, F_2)

F_3 = HU.fayRiddelCase3(h_dOverH_0e)
qDot_w3 = HU.fayRiddellHeatFlux(k, Pr_w, rho_0e, mu_0e, rho_w, mu_w, dudx_e, h_0e, h_w, F_3)