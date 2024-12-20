clear; close all; clc; addpath('..')
HU = HypersonicsUtilities;

%% 1
R_N = 0.55; % m
B = 6000; % Pa
c_D = 1.4;
k = 2e-4/sqrt(R_N)
LDRatio = 1.4;
v_0 = 7950; % m/s
h_0 = 100e3; % m

%% 1a
qDot_max = HU.maxHeatingLifting(k, B, LDRatio, v_0)
v_qmax = HU.velocityAtMaxHeatingLifting(v_0)
rho_qmax = (qDot_max/(k*v_qmax^3))^2

rho_sl = HU.airDensity_seaLevel;
aTil = HU.inverseScaleHeight;
h = log(rho_qmax/rho_sl)/-aTil

%% 1b
qDot_Detra = HU.detraEtAllHeatFlux(rho_qmax, v_qmax, R_N)

%% 1c
sigma = HU.stefanBoltzmannConstant;
e = 0.85;
T_w = (qDot_Detra/(e*sigma))^(1/4)

T_wLimit = 2000; % K
qDot_limit = HU.radiationHeatFlux(T_wLimit, e)
saftyFactor = qDot_limit/qDot_Detra

%% 1ei
h = 70e3; % m
L = 40; % m
T_inf = 219.585; % K
rhoRatio = 6.7616e-5;
rho_inf = rho_sl*rhoRatio
mu_inf = HU.sutherlandsLawViscosity(T_inf)
Re_inf = HU.reynoldsNumber(rho_inf, v_qmax, L, mu_inf)

%% 1fi
q = .5*rho_inf*v_qmax^2
a_inf = HU.speedOfSoundAir(T_inf)
M_inf = v_qmax/a_inf

%% 1fii
atm = HU.pascalsPerAtmosphere;
phi = q/atm