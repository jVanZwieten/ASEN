clear;clc;close all;
addpath("..\..")

c   = 299792458;    % SI speed of light, m/s
omegaE = 7.2921151467e-5;
f_L1 = 1575.42e6;
f_L2 = 1227.6e6;

rVecEcef_nist = [-1288398.574; -4721696.936; 4078625.349];
load('data_rinex_NIST00USA_R_2024_239.mat')
navfilename='brdc2390.24n';
[gps_ephem,ionoparams]=read_GPSbroadcast(navfilename);

%% 1a
initialEpoch = datetime(2024, 8, 26, 0, 0, 0);
rinexTt_firstEpoch = data_rinex.GPS(data_rinex.GPS.Time == initialEpoch, :);

[epochInitial_gpsWeek, epochInitial_timeOfWeek, epochInitial_gpsWeekMod1024] = cal2gps(initialEpoch);

varNames = {'prn', 'rVecEcef', 'vVecEcef'};
varTypes = {'double', 'cell', 'cell'};
satStateT = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);

varNames = {'prn', 'residual', 'f_doppler', 'azimuth', 'elevation'};
varTypes = {'double', 'double', 'double', 'double', 'double'};
satDataT = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);

for i = 1:size(rinexTt_firstEpoch, 1)
    prn = rinexTt_firstEpoch.SatelliteID(i);
    [health, rVecEcef_sat, vVecEcef_sat, clockCorrection, reletivityCorrection, tgd0] = eph2pvt(gps_ephem, [epochInitial_gpsWeek, epochInitial_timeOfWeek], prn);
    rVecEcef_sat = rVecEcef_sat(:);
    vVecEcef_sat = vVecEcef_sat(:);

    [Az, El, Range] = compute_azelrange(rVecEcef_nist', rVecEcef_sat');

    C1W = rinexTt_firstEpoch.C1W(i);
    C2W = rinexTt_firstEpoch.C2W(i);

    delay_iono = ionosphericDelayL1(C1W, C2W);
    delay_tropo = troposphereDelay(initialEpoch, rVecEcef_nist', El);
    residual = C1W - Range - clockCorrection - reletivityCorrection - delay_iono - delay_tropo;

    f_doppler = dopplerFrequency(rVecEcef_sat', vVecEcef_sat', rVecEcef_nist, f_L1);

    satSate = [satStateT; {prn, rVecEcef_sat, vVecEcef_sat}];
    satDataT = [satDataT; {prn, residual, f_doppler, Az, El}];
end

%% 1b


function delay_ionoL1 = ionosphericDelayL1(C1C, C2W)
    f_L1 = 1575.42e6;
    f_L2 = 1227.6e6;
    totalElectronCount = f_L1^2*f_L2^2/40.3/(f_L1^2-f_L2^2)*(C2W-C1C);
    delay_ionoL1 = 40.3*totalElectronCount/f_L1^2;
end

function delay_tropo = troposphereDelay(epoch, rVecEcef_rx, elevation)
    epoch = Utilities.dayOfYearTimeFraction(epoch);
    rLla_rx = ecef2lla(rVecEcef_rx);
    orth_height = rLla_rx(3)-geoidheight(rLla_rx(1), rLla_rx(2));
    [delay_tropo] = UNB3M(deg2rad(rLla_rx(1)), orth_height, epoch, deg2rad(elevation));
end

function f_doppler = dopplerFrequency(rMEcef_tx, vMEcef_tx, rVEcef_rx, f_base)
    c = Utilities.speedOfLight;

    rHatEcef_txRx = Utilities.unitVector(rVEcef_rx - rMEcef_tx');
    f_doppler = (dot(vMEcef_tx', rHatEcef_txRx)*f_base/c)';
end