format shortG

HU = HypersonicsUtilities;
rho_sl = HU.airDensity_seaLevel;
p_sl = HU.airPressure_seaLevel;
gamma_air = HU.ratioSpecificHeats_air;
R_air = HU.gasConstant_air;

%% 3
L = 32.77; % m
u_1 = 5100; % m/s
h = 60; % km

% 3a
T_1 = 253.114; % K
rho_1 = 2.5280e-4*rho_sl
p_1 = 2.1671e-4*p_sl
a_1 = 315.07; % m/s

M_1 = u_1/a_1
mu_1 = HU.sutherlandsLawViscosity(T_1)
Re_1 = HU.reynoldsNumber(rho_1, u_1, L, mu_1)

% 3b
rhoURatio_12 = HU.densityVelocity12Ratio(gamma_air, M_1)
rho_2 = rhoURatio_12*rho_1
u_2 = u_1/rhoURatio_12
pRatio_12 = HU.pressureRatio(gamma_air, M_1)
p_2 = p_1*pRatio_12
M_2 = HU.machAfterShock(gamma_air, M_1)
T_2 = p_2/(rho_2*R_air)

a_2 = HU.speedOfSound(gamma_air, R_air, T_2)
u_2check = M_2*a_2

%% 3c
p_2 = p_1 + rho_1*u_1^2
p_2/HU.pascalsPerAtmosphere
c_p = HU.specificHeatAtConstantPressure(gamma_air, R_air)
h_1 = c_p*T_1
h_2 = HU.enthalpyTotal(h_1, u_1)

T_2 = 5200; % K
z_2 = 1.3;
rho_2 = HU.densityRealGasEffect(p_2, z_2, R_air, T_2)
u_2 = rho_1*u_1/rho_2

p_2 = p_2 - rho_2*u_2^2
p_2/HU.pascalsPerAtmosphere
h_2 = h_2 - .5*u_2^2

z_2 = 1.298;
rho_2 = HU.densityRealGasEffect(p_2, z_2, R_air, T_2)
u_2 = rho_1*u_1/rho_2

p_23 = p_1 + rho_1*u_1^2 - rho_2*u_2^2
(p_2 - p_23)/p_2

% p_ex1 = 5.1526e-05*p_sl
% rho_ex1 = 6.7616e-05*rho_sl
% p_ex2 = p_ex1 + rho_ex1*7000^2
% p_ex2/HU.pascalsPerAtmosphere

%% 3d
p_23atm = p_23/HU.pascalsPerAtmosphere
mu_2s = HU.sutherlandsLawViscosity(T_2)
mu_2 = 1.1*mu_2s
k_2s = HU.sutherlandsLawThermalConductivity(T_2)
k_2 = 12*k_2s

%% 3e
p_02 = p_2 + rho_2*u_2^2
p_02 = p_02/HU.pascalsPerAtmosphere
a = 221;
b = .029;
t_v = exp(a*(T_2^(-1/3) - b) - 18.42)/p_02
x_v = u_2*t_v

%% 3f
rho_dN2 = 130000;
theta_dN2 = 113000;
aSqrdOver1MinusA = rho_dN2/rho_2*exp(-theta_dN2/T_2)
alpha = roots([1 aSqrdOver1MinusA -aSqrdOver1MinusA])
a = 1;
b = aSqrdOver1MinusA;
c = -aSqrdOver1MinusA;
alpha = (-b + sqrt(b^2 - 4*a*c))/(2*a)
alpha = (-b - sqrt(b^2 - 4*a*c))/(2*a)