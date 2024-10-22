clear
close all
addpath('..')

%% 1
mu_earth = 398600.5; % km^3/s^2
c = 3e8; % m/s

r_gps = 26560; % km
h_leo = 800; % km
i = 55; % degrees

r_earth = 6380; % km
r_leo = r_earth + h_leo

f_gps = 1575.42e6; % Hz
f_leo = 12e9; % Hz

%% 1a
v_gps = sqrt(mu_earth/r_gps)
v_gps = v_gps*1000;
v_leo = sqrt(mu_earth/r_leo)
v_leo = v_leo*1000;

%% 1b
h_gps = r_gps - r_earth

%% 1c
r_maxGps = sqrt(r_gps^2 - r_earth^2)
r_maxLeo = sqrt(r_leo^2 - r_earth^2)

%% 1e
vRhat_gps = v_gps*r_earth/r_gps
f_dopplerGps = f_gps*vRhat_gps/c

vRhat_leo = v_leo*r_earth/r_leo
f_dopplerLeo = f_leo*v_leo/c

%% 1g
r_zLeo = r_leo*sind(55)

%% 2
rVLla_tx = [deg2rad(40); deg2rad(30); 100]; % rad rad m
f_tx = 1575.42e6; % Hz

%% 2a
format long
r_tx = r_earth*1000 + rVLla_tx(3)
rZ_tx = r_tx*sin(rVLla_tx(1))
rXY_tx = r_tx*cos(rVLla_tx(1))
rX_tx = rXY_tx*cos(rVLla_tx(2))
rY_tx = rXY_tx*sin(rVLla_tx(2))

rVEcef_tx = [rX_tx; rY_tx; rZ_tx]
assert((norm(rVEcef_tx) - r_tx) < .01)

%% 2b
format shortG
lat = rVLla_tx(1);
long = rVLla_tx(2);
C_Ecef2Enu = [-sin(long) cos(long) 0
-sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat)
cos(lat)*cos(long) cos(lat)*sin(long) sin(lat)]

%% 2c
el_rxTx = atand(1/9)
rg_rxTx = norm([100 900])

%% 2e
format longG
rVEnu_rxTx = [900; 0; 100]; % m
rVEcef_rx = rVEcef_tx - C_Ecef2Enu'*rVEnu_rxTx
norm(rVEcef_tx - rVEcef_rx)

%% 2f
vVEnu_B = [30; 0; 0]; % m
rHatEnu_txRx = -rVEnu_rxTx/norm(rVEnu_rxTx)
rDot_bTx = dot(vVEnu_B, rHatEnu_txRx)

f_doppler = f_tx*rDot_bTx/c

%% 2g
vVEnu_B2 = [0; 30; 0]; % m/s
rDot_bTx2 = dot(vVEnu_B2, rHatEnu_txRx)

%% 3b
registerA = ones(1, 4);
codeA = codeFromRegister(registerA, [4 3], 15)
codeB = codeFromRegister(registerA, [4 1], 15)

%% 3c
R_A = correlate(codeA, codeA)
R_B = correlate(codeA, codeB)

figure
plot(1:15, R_B)
xlabel("n")
ylabel("correlation R(n)")
title("Correlation of Codes A & B")

%% 3e
chipRate = 3e6; % chips/s
chipTime = 1/chipRate
chipLength = c/chipRate

function code = codeFromRegister(register, keys, iterations)
    code = zeros(1, iterations);
    for i = 1:iterations
        code(i) = registerOutput(register(end));
        register = [xor(register(keys(1)), register(keys(2))) register(1:end - 1)];
    end
end

function output = registerOutput(registerInput)
    output = registerInput*-2 + 1;
end

function R = correlate(codeA, codeB)
    assert(length(codeA) == length(codeB))
    R = NaN(1, length(codeA));
    for n = 1:length(R)
        codeB_n = Utilities.leftShift(codeB, n);
        R(n) = sum(codeA.*codeB_n);
    end
end