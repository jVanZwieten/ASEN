addpath("..")
clear

HU = HypersonicsUtilities;

%% 1
R_N = 4.8; % m
R_B = 2.25; % m
A = pi*R_B^2;
m = 7500; % kg
c_D = 1.6;
L = 0;
T = 0;
k = 2e-4/sqrt(R_N)
v_0 = 11100; % m/s
h_0 = 100e3; % m
gamma_0 = deg2rad(-7);

%% 1a
beta = -2*m*HU.inverseScaleHeight*sin(gamma_0)/(1.6*A)

qDot_max = HU.maxHeatingBallistic(k, beta, v_0)
g_max = HU.maxG(v_0, gamma_0)

%% 2b
T_1 = 264.164; % K
P_1 = 149.1; % Pa
rho_1 = 1.9662e-3; % kg/m^3
u_1 = 5000; % m/s
h = 45e3; % m

P_2 = P_1 + rho_1*u_1^2
c_p = 1010;
h_2 = c_p*T_1 + .5*u_1^2
P_2atm = P_2/HU.pascalsPerAtmosphere
h_2/10^6

z = 1.3;
T_2 = 5400;

rho_2 = HU.densityRealGasEffect(P_2, z, HU.gasConstant_air, T_2)
u_2 = rho_1*u_1/rho_2

%% 2c
mu = HU.sutherlandsLawViscosity(T_2)

k = HU.sutherlandsLawThermalConductivity(T_2)
k = 7.5*k