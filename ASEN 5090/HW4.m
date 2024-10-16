if ~exist("sp3", "var")
    sp3 = read_sp3('COD0MGXFIN_20242390000_01D_05M_ORB.SP3');
end
if ~exist("rinexGPS", "var")
    rinexGPS = load('data_rinex_NIST00USA_R_2024_239.mat').data_rinex.GPS;
end
if ~exist("ephemerisData", "var") || ~exist("ionoParams", "var")
    [ephemerisData, ionoParams] = read_GPSbroadcast('brdc2390.24n');
end

clearvars -except sp3 rinexGPS ephemerisData ionoParams
close all
addpath('..')

weekOfYear =  2329;
weekSeconds_0 = 86400;
dayOfYear = 239;
dayStart = datetime(2024, 8, 26, 0, 0, 0, "TimeZone", "UTC");
dayEnd = datetime(2024, 8, 26, 24, 0, 0, "TimeZone", "UTC");

GU = GnssUtilities;

%% 1
rinexTt_prn2_16 = rinexGPS((rinexGPS.SatelliteID == 2 | rinexGPS.SatelliteID == 16) & ~isnan(rinexGPS.C1C) & ~isnan(rinexGPS.L1C) & ~isnan(rinexGPS.C2W) & ~isnan(rinexGPS.L2W), {'SatelliteID' 'C1C' 'L1C' 'C2W' 'L2W'});
delayV_ionosphericL1Gps = GU.ionosphericDelayL1L2(GU.wavelength_L1*rinexTt_prn2_16.L1C, GU.wavelength_L2*rinexTt_prn2_16.L2W);
check = GU.ionosphericDelayL1L2(rinexTt_prn2_16.C1C, rinexTt_prn2_16.C2W);

rinexTt_prn2_16.('ionosphericDelay') = delayV_ionosphericL1Gps;
rinexTt_prn2_16.('multipathError') = rinexTt_prn2_16.C1C - GU.wavelength_L1*rinexTt_prn2_16.L1C + 2*rinexTt_prn2_16.ionosphericDelay;

multipathErrorTt_2 = rinexTt_prn2_16(rinexTt_prn2_16.SatelliteID == 2, {'multipathError'});
multipathErrorTt_16 = rinexTt_prn2_16(rinexTt_prn2_16.SatelliteID == 16, {'multipathError'});

figure
hold on
plot(multipathErrorTt_2.Time, multipathErrorTt_2.multipathError)
plot(multipathErrorTt_16.Time, multipathErrorTt_16.multipathError)

%% 2
rVEcef_nist = [-1288398.574; -4721696.936; 4078625.349];
rVGeod_nist = GU.Ecef2Geodetic(rVEcef_nist);

tt_prn8 = rinexGPS(rinexGPS.SatelliteID == 8 & all(~isnan([rinexGPS.C1C rinexGPS.C2W]), 2), {'C1C' 'C2W'});
[azElRgM_trx, rMEcef_prn8, vMEcef_prn8, clockCorrectionV_prn8, relativisticCorrectionV_prn8] = GU.gpsStateXmitOrigin(rVEcef_nist, ephemerisData, 8, tt_prn8.Time);
tt_prn8.('ExpectedRange') = azElRgM_trx(3, :)';
tt_prn8.('ClockBias') = clockCorrectionV_prn8;
tt_prn8.('RelativisticEffect') = relativisticCorrectionV_prn8;
tt_prn8.('IonosphereDelay') = GU.ionosphericDelayL1L2(tt_prn8.C1C, tt_prn8.C2W);

tropoDelay = zeros(size(tt_prn8, 1), 1);
for i = 1:length(tropoDelay)
    tropoDelay(i) = UNB3M(rVGeod_nist(1), rVGeod_nist(3), dayOfYear, azElRgM_trx(2, i));
end
tt_prn8.('TroposphereDelay') = tropoDelay;

tt_prn8.('Residual') = tt_prn8.C1C - tt_prn8.ExpectedRange - tt_prn8.ClockBias - tt_prn8.RelativisticEffect - tt_prn8.IonosphereDelay - tt_prn8.TroposphereDelay;

tt_prn8 = Utilities.breakTimeTable(tt_prn8, datetime(2024, 08, 26, 10, 0, 0));

figure
stackedplot(tt_prn8)

%% 3
dopplerTt_prn8 = rinexGPS(rinexGPS.SatelliteID == 8 & all(~isnan([rinexGPS.C1C rinexGPS.C2W]), 2), {'D1C'});
dopplerTt_prn8.('DopplerShift') = GU.dopplerFrequency(rMEcef_prn8, vMEcef_prn8, rVEcef_nist, GU.frequency_L1);
dopplerTt_prn8.('delta') = (dopplerTt_prn8.DopplerShift - dopplerTt_prn8.D1C)./dopplerTt_prn8.DopplerShift;
dopplerTt_prn8 = Utilities.breakTimeTable(dopplerTt_prn8, datetime(2024, 08, 26, 10, 0, 0));
figure
stackedplot(dopplerTt_prn8)