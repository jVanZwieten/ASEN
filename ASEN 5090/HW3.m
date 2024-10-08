clear
close all
addpath("..")

weekNumber =  2329;
weekSeconds_0 = 86400;

%% 1
sp3 = read_sp3('COD0MGXFIN_20242390000_01D_05M_ORB.SP3');
%   1 - Week number
%   2 - TOW (s)
%   3 - PRN
%   4-6 - XYZ Position (km)
%   7 - Clock bias (microsec)
REcf_Gps =  sp3(sp3(:,8) == 1, 2:6)';
REcf_Gps(1, :) = (REcf_Gps(1, :) - weekSeconds_0) / (60*60);
REcf_Gps(3:5, :) = 1000*REcf_Gps(3:5, :);
% 1 - TOW (h)
% 2 - PRN
% 3-5 - ECEF XYZ (m)

REcef_Nist = [-1288398.574; -4721696.936; 4078625.349];
minVisibleElevation = deg2rad(5);
SatellitesInView = GnssUtilities.SatellitesInViewAtEcef(REcf_Gps(1:5, :), REcef_Nist, minVisibleElevation);

figure
plot(SatellitesInView(1, :), SatellitesInView(2, :), 'x');
xlim([0 24])
xticks(0:2:24)
xlabel("Time (h)")
ylim([1 32])
yticks(1:32)
ylabel("GPS PRN")
title("GPS Satellite Visibility through 26AUG24")
grid on

rinexGPS = load('data_rinex_NIST00USA_R_2024_239.mat').data_rinex.GPS;
prn10 = rinexGPS(rinexGPS.SatelliteID == 10, {'C1C', 'C1W'});

figure
stackedplot(prn10, "DisplayLabels", ["C1C (m)", "C1W (m)"])
title("PRN10 Psedoranges through Transit")

%% 2
[~, ionoParams] = read_GPSbroadcast('brdc2390.24n');
alphas = ionoParams(1:4);
betas = ionoParams(5:8);

% test
% RGeod_Nist = GnssUtilities.Ecef2Geodetic(REcef_Nist);
% t_0 = 0;
% t_f = Utilities.days2s(1);
% long_IPP = GnssUtilities.l
% localTime_0 = GnssUtilities.localTime(t_0, RGeod_Nist(2))
% localTime_final = GnssUtilities.localTime(t_f, RGeod_Nist(2))

% latGeomag_IPP0 = GnssUtilities.latitudeGeomagnetic()