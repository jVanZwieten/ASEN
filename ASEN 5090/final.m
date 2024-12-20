clc; clear; close all; addpath('..')
GU = GnssUtilities;
c = 3e8; % m/s

%% 1
prn = [1:5]';
x = [-8 -6 0 1 3]';
y = [0 6 3 1 7]';
pseudorange = [17 15.5 15.9 17.4 20.2]';
satelliteData = table(prn, x, y, pseudorange);

%% 1a
XVec_guess = [mean(satelliteData.x); mean(satelliteData.y); 0];
XVec_RxEstimate = XVec_guess;

%% 1bc
positionMat_Tx = [satelliteData.x satelliteData.y];
dX_i = inf;

while(norm(dX_i) > 1e-4)
    positionVec_RxEstimate = XVec_RxEstimate(1:2);
    pseudorange_TxRxEstimate = range(positionVec_RxEstimate, positionMat_Tx) + XVec_RxEstimate(3);
    dPseudorange = satelliteData.pseudorange - pseudorange_TxRxEstimate;
    G = geometryMatrix(positionVec_RxEstimate, positionMat_Tx);
    
    dX_i = G\dPseudorange;
    XVec_RxEstimate = XVec_RxEstimate + dX_i;
end
XVec_RxEstimate

%% 1e
H = inv(G'*G);
HDOP = sqrt(H(1,1) + H(2,2))

%% 1f
satelliteDataSubset = satelliteData([1 4 5], :);
XVec_RxEstimate = XVec_guess;
positionMat_Tx = [satelliteDataSubset.x satelliteDataSubset.y];

dX_i = inf;
while(norm(dX_i) > 1e-4)
    positionVec_RxEstimate = XVec_RxEstimate(1:2);
    pseudorange_TxRxEstimate = range(positionVec_RxEstimate, positionMat_Tx) + XVec_RxEstimate(3);
    dPseudorange = satelliteDataSubset.pseudorange - pseudorange_TxRxEstimate;
    G = geometryMatrix(positionVec_RxEstimate, positionMat_Tx);
    
    dX_i = G\dPseudorange;
    XVec_RxEstimate = XVec_RxEstimate + dX_i;
end
XVec_RxEstimate
H = inv(G'*G);
HDOP = sqrt(H(1,1) + H(2,2))

%% 2
alt = 1000e3; % m
R_earth = 6378e3; % m

%% 2a
f_c = 1190e6; % Hz
wavelength_c = c/f_c % m

%% 2b
Ls = GU.spaceLoss(alt, wavelength_c) % dB

f_L1 = 1575.42e6; % Hz
wavelength_L1 = c/f_L1 % m
alt_GPS = 20000e3; % m
Ls_GPS = GnssUtilities.spaceLoss(alt_GPS, wavelength_L1) % dB

%% 2c
R_orbit = alt+R_earth;
el = 0;
alpha = asin(R_earth*sin(el + pi/2)/R_orbit)
rad2deg(alpha)

offNadirMultiplier=0.75;
offNadirElement = 10*log10(offNadirMultiplier);

R_horizon = R_orbit*cos(alpha)
Ls_horizon = GnssUtilities.spaceLoss(R_horizon, wavelength_c) % dB
Ls_total = -Ls_horizon + offNadirElement

%% 2d
TEC = 10e16;
I_rho = GU.ionosphericDelay(TEC, f_c) % m

function G = geometryMatrix(XVec, XMat_sats)
    G = XVec' - XMat_sats;
    R = range(XVec, XMat_sats);
    G = G ./ R;
    G = [G ones(size(G, 1), 1)];
end

function r = range(XVec, XMat)
    r = sqrt(sum((XVec' - XMat).^2, 2));
end