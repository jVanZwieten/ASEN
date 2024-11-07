addpath("..")
clear

HU = HypersonicsUtilities;

%% 2a
T = 3500; % K
p = 1; % atm
mu = HU.sutherlandsLawViscosity(T)
k = HU.sutherlandsLawThermalConductivity(T)
k = 5*k

%% 2d
rho = .1;
alphaOverAlpha1 = HU.massFractionFunction(rho, 150000, 59500, T)
alpha = roots([1 alphaOverAlpha1 -alphaOverAlpha1])
alpha = alpha(2)
C_O = alpha*.2
C_O2 = (1 - alpha)*.2

%% 2e
R_N2 = HU.speciesGasConstant(28.0134)
R_O2 = HU.speciesGasConstant(31.9988)
R_O = HU.speciesGasConstant(15.9994)

%% 2f
c_pN2 = HU.specificHeatConstantPressureDiatomic(R_N2, 3390, T)
c_pN2/R_N2
c_pO2 = HU.specificHeatConstantPressureDiatomic(R_O2, 2270, T)
c_pO2/R_O2
c_pO = HU.specificHeatConstantPressureMonoatomic(R_O)
c_pO/R_O

%% 2g
c_pMix = sum([.8*c_pN2 C_O2*c_pO2 C_O*c_pO])

%% 3
v_0 = 5000; % m/s
h_0 = 20e3; % m
beta = .5; % kg/m^3
A = .1*.1; % m^2
phi = deg2rad(30)
T = 250; % K

rho_1 = HU.airDensity(20e3)
v_1 = v_0*HU.velocityRatio(rho, beta)
gamma = 1.4;
a_1 = HU.speedOfSound(gamma, HU.gasConstant_air, T);
M_1 = v_1/a_1
rhoRatio = HU.densityVelocity12Ratio(gamma, M_1)
rho_2 = rho_1*rhoRatio
v_2 = v_1/rhoRatio
M_2 = HU.machAfterShock(gamma, M_1)

%% 2a
R_nose = 1; %m
k = 2e-7/sqrt(R_nose)
LtoD = 1.2;
B = 2500;
v_0 = 7950;
h_0 = 100e3;
q_peak = HU.maxHeatingLifting(k, B, LtoD, v_0)

%% 2b
v = HU.velocityAtMaxHeatingLifting(7950)
rho = (1 - v^2/v_0^2)/(LtoD*v^2/2/B)
h = HU.altitudeFromAirDensity(rho)

%% 4
u_inf = 6000; % m/s
h = 40e3; % m
L = 5; %m
gamma = 1.4;
%% 4a
T_inf = 250.35; % K
rho_inf = 3.2618e-3*HU.airDensity_seaLevel
p = rho_inf*HU.gasConstant_air*T_inf
a = HU.speedOfSound(gamma, HU.gasConstant_air, T)
M = u_inf/a
mu = HU.sutherlandsLawViscosity(T_inf)
Re = HU.reynoldsNumber(rho_inf, u_inf, L, mu)
%% 4b
rhoRatio = HU.densityVelocity12Ratio(gamma, M)
rho_2 = rho_inf*rhoRatio
u_2 = u_inf/rhoRatio
p_2 = p*HU.pressureRatio(gamma, M)
M_2 = HU.machAfterShock(gamma, M)
T_2 = p_2/rho_2/HU.gasConstant_air
%% 4c
p_2 = p + rho_inf*u_inf^2
p_2 = p_2/HU.pascalsPerAtmosphere
c_p = 1010;
h_2 = c_p*T_inf + .5*u_inf^2
z = 1.39;
T_2 = 6500;
rho_2 = p_2*HU.pascalsPerAtmosphere/z/HU.gasConstant_air/T_2
u_2 = rho_inf*u_inf/rho_2
%% 4d
mu = HU.sutherlandsLawViscosity(T_2)
mu = 1.1*mu
Re_2 = HU.reynoldsNumber(rho_2, u_2, L, mu)