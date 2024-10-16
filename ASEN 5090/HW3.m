close all
addpath("..")

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

weekOfYear =  2329;
weekSeconds_0 = 86400;
dayOfYear = 239;
dayStart = datetime(2024, 8, 26, 0, 0, 0, "TimeZone", "UTC");
dayEnd = datetime(2024, 8, 26, 24, 0, 0, "TimeZone", "UTC");

dayHourLimit = [0 24];
dayHourLimitDatetime = [dayStart dayEnd];

%% 1
%   1 - Week number
%   2 - TOW (s)
%   3 - PRN
%   4-6 - XYZ Position (km)
%   7 - Clock bias (microsec)
REcf_Gps =  sp3(sp3(:,8) == 1, 2:6)';
REcf_Gps(1, :) = (REcf_Gps(1, :) - weekSeconds_0);
REcf_Gps(3:5, :) = Utilities.km2m(REcf_Gps(3:5, :));
% 1 - TOW (s)
% 2 - PRN
% 3-5 - ECEF XYZ (m)

REcef_Nist = [-1288398.574; -4721696.936; 4078625.349];
minVisibleElevation = deg2rad(5);
tAzElRg_Gps = [REcf_Gps(1:2, :); GnssUtilities.AzElRangeToObserverEcef(REcef_Nist, REcf_Gps(3:5, :))];
SatellitesInView = GnssUtilities.SatellitesInView(tAzElRg_Gps, minVisibleElevation);

figure
plot(Utilities.s2h(SatellitesInView(1, :)), SatellitesInView(2, :), 'x');
xlim([0 24])
xticks(0:2:24)
xlabel("Time (h)")
ylim([1 32])
yticks(1:32)
ylabel("GPS PRN")
title("GPS Satellite Visibility through 26AUG24")
grid on

rhos_prn10 = rinexGPS(rinexGPS.SatelliteID == 10, {'C1C', 'C1W'});

figure
stackedplot(rhos_prn10, "DisplayLabels", ["C1C (m)", "C1W (m)"])
title("PRN10 Psedoranges through Transit")

%% 2
alphas = ionoParams(1:4);
betas = ionoParams(5:8);

RGeod_Nist = GnssUtilities.Ecef2Geodetic(REcef_Nist);
altitude_Nist = 1615;
tAzElRg_prn8 = SatellitesInView([1 3:5], SatellitesInView(2, :) == 8);

Klobuchar = zeros(2, size(tAzElRg_prn8, 2));
for i = 1:size(Klobuchar,  2)
    AzEl = tAzElRg_prn8(:, i);
    Klobuchar(:, i) = [AzEl(1); GnssUtilities.Klobuchar(RGeod_Nist, AzEl(2:3), AzEl(1), alphas, betas)*Utilities.speedOfLight];
end

KlobucharPlotData = [Utilities.s2h(Klobuchar(1, :)); Klobuchar(2, :)];
KlobucharPlotData = [Utilities.breakSeries(KlobucharPlotData, 10)]; % separates 2 'islands' of data without connecting trendline

figure
plot(KlobucharPlotData(1, :), KlobucharPlotData(2, :))
xlabel("t_{GPS} (h)")
xlim(dayHourLimit)
ylabel("Delay_{Ionosphere} (m)")
title("Ionospheric Delay for PRN 8 L1 over NIST, Klobuchar Model")

rhos_prn8 = rinexGPS(rinexGPS.SatelliteID == 8 & ~isnan(rinexGPS.C1C) & ~isnan(rinexGPS.C2L), {'C1C' 'C2L'});
I_prn8 = GnssUtilities.ionosphericDelayL1L2(rhos_prn8.C1C, rhos_prn8.C2L);

empiricalIonosphereData = [Utilities.hourOfDay(rhos_prn8.Time)'; I_prn8'];
empiricalIonosphereData = Utilities.breakSeries(empiricalIonosphereData, 10);

figure
plot(empiricalIonosphereData(1, :), empiricalIonosphereData(2, :))
xlabel("t_{GPS} (h)")
xlim(dayHourLimit)
ylabel("Delay_{Ionosphere} (m)")
title("Ionospheric Delay for PRN 8 L1 over NIST, Empirical")

%% 3
reasonableZenithTopoDelay_NIST = 2; % m
topoDelaySimple_prn8 = GnssUtilities.toposphericDelaySimple(reasonableZenithTopoDelay_NIST, tAzElRg_prn8(3, :));
topoDelayPlotData = [Utilities.s2h(tAzElRg_prn8(1, :)); topoDelaySimple_prn8];
topoDelayPlotData = Utilities.breakSeries(topoDelayPlotData, 10);

topoXLimit = [1 21];
topoYLimit = [0 22];
figure
plot(topoDelayPlotData(1, :), topoDelayPlotData(2, :))
xlabel("Time (h)")
xlim(topoXLimit)
ylabel("Delay (m)")
ylim(topoYLimit)
title("Tropospheric Delay, PRN8, Simple Model")

unb3mDelay = zeros(3, size(tAzElRg_prn8, 2));
unb3mDelay(1, :) = Utilities.s2h(tAzElRg_prn8(1, :));
for i = 1:size(tAzElRg_prn8, 2)
    unb3mDelay(2, i) = UNB3M(RGeod_Nist(1), altitude_Nist, dayOfYear, tAzElRg_prn8(3, i));
end
unb3mDelay = Utilities.breakSeries(unb3mDelay, 10);

figure
plot(unb3mDelay(1, :), unb3mDelay(2, :))
xlabel("Time (h)")
xlim(topoXLimit)
ylabel("Delay (m)")
ylim(topoYLimit)
title("Tropospheric Delay, PRN8, UNB3M model")

%% 4
T = [weekOfYear*ones(size(tAzElRg_prn8, 2), 1), tAzElRg_prn8(1, :)'+weekSeconds_0];

[~, REcef_prn8] = eph2pvt(ephemerisData, T, 8);
AzElFromEphemeris_prn8 = GnssUtilities.AzElRangeToObserverEcef(REcef_Nist, REcef_prn8');

T = GnssUtilities.gps2datetime(T);
AzElTransitTimeCorrected_prn8 = GnssUtilities.AzElRangeToObserverEcefAtTimeReception(REcef_Nist, ephemerisData, 8, T);

rangett = timetable(T);
rangett.("Range_uncorrected") = AzElFromEphemeris_prn8(3, :)';
rangett.("Range_corrected") = AzElTransitTimeCorrected_prn8(3, :)';
rangett = Utilities.appendTimetableRow(rangett, datetime(2024, 8, 26, 10, 0, 0, 'Timezone', 'UTC'), {NaN NaN});

figure
stackedplot(rangett)
xlim(dayHourLimitDatetime)
title("PRN8 Range to NIST at Transmission Time")

figure
plot(rangett.T, rangett.Range_corrected - rangett.Range_uncorrected);
xlabel("Time")
xlim(dayHourLimitDatetime)
ylabel("Range Correction (m)")
title("Range Correction for Transit Time, PRN8")

testDateTimes = [dayStart; dayEnd];
[~, REcef_test] = eph2pvt(ephemerisData, GnssUtilities.datetime2Gps(testDateTimes), 8);
AzElRg_test = GnssUtilities.AzElRangeToObserverEcef(REcef_Nist, REcef_test');
AzElRg_testCorrected = GnssUtilities.AzElRangeToObserverEcefAtTimeReception(REcef_Nist, ephemerisData, 8, testDateTimes);
testRanges = timetable(testDateTimes, AzElRg_test(3, :)', AzElRg_testCorrected(3, :)');