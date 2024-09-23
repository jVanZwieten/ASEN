HU = HypersonicsUtilities;

%% 2
R_A = .05; % m
R_B = .9; % m
A = pi*R_B^2; % m^2
m = 3400; % kg
c_D = .4; % drag coefficient
k = 2e-4/sqrt(R_B)
v_0 = 4500; % m/s
h_0 = 100; % km
gamma = deg2rad(-70); % flight path angle, degrees

%% 2a
beta = HU.beta(m, gamma, c_D, A)
heating_max = HU.maxHeating(k, beta, v_0)
v_heatingMax = HU.velocityRatio_heatingMax*v_0
rho_heatingMax = (heating_max/k/v_heatingMax^3)^2
h_heatingMax = HU.altitude(rho_heatingMax)

decel_max = HU.maxG(v_0, gamma)
rho_gMax = beta/2
h_gMax = HU.altitude(rho_gMax)
v_gMax = HU.velocityRatio_decelMax*v_0