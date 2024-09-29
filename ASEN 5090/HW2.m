addpath('..')
clear
close all

%% 1-1
c = 2.99792458e8; % m/s

gps_almanac = read_GPSyuma('almanac.yuma.week0281.147456.txt', 2);
gps_ephem = read_GPSbroadcast('brdc2390.24n');
sp3 = read_sp3('COD0MGXFIN_20242390000_01D_05M_ORB.SP3'); % considered to be truth

%% 1-2
prn = 19;

Orbit_prn19true =  sp3(sp3(:,3) == prn & sp3(:,8) == 1, [1:2 4:7]);
%   1 - Week number
%   2 - TOW (s)
%   4-6 - XYZ Position (km)
%   7 - Clock bias (microsec)
Orbit_prn19true = Orbit_prn19true(1: end - 1, :); % omit erroneous final epoch
Orbit_prn19true(:, 3:5) = Orbit_prn19true(:, 3:5)*1000; % km -> m
Orbit_prn19true(:, 6) = Orbit_prn19true(:, 6)/1e6*c; % ms -> m

T = Orbit_prn19true(:, 1:2);
[~, Orbit_prn19Almanac, ClockBias_prn19Almanac] = alm2pos(gps_almanac, T, prn);
dAlmanac = [T  ([Orbit_prn19Almanac ClockBias_prn19Almanac] - Orbit_prn19true(:, 3:end))];
plotGPSDifference(dAlmanac(:, 2:end)', "Almanac Prediction - SP3, PRN 19")

[~, Orbit_prn19Ephemeris, ~, ClockBias_prn19Ephemeris] = eph2pvt(gps_ephem, T, prn);
dEphemeris = [T ([Orbit_prn19Ephemeris ClockBias_prn19Ephemeris] - Orbit_prn19true(:, 3:end))];
plotGPSDifference(dEphemeris(:, 2:end)', "EphemerisPrediction - SP3, PRN 19")

%% 1-2
[~, Orbit_prn19EphemerisExtrapolate, ~, ClockBias_prn19EphemerisExtrapolate] = eph2pvtFromFirst(gps_ephem, T, prn);
dEphemerisExtrapolate = [T ([Orbit_prn19EphemerisExtrapolate ClockBias_prn19EphemerisExtrapolate] - Orbit_prn19true(:, 3:end))];
plotGPSDifference(dEphemerisExtrapolate(:, 2:end)', "EphemerisPrediction(Extrapolated) - SP3, PRN 19")

%% 2-1
R_observerBoulder = GnssUtilities.Geodetic2EcefSurface(deg2rad(40.010228), deg2rad(-105.243872));
noonMdt = 151200; % week seconds
t = [2329 noonMdt]; % GPS week, week second

GpsInView = sp3(sp3(:, 8) == 1 & sp3(:, 2) == noonMdt, 3:6)';
GpsInView = [GpsInView(1, :); GpsInView(2:4, :)*1000];
for i = 1:size(GpsInView, 2)
    [az, el, range] = computeAzElRange(R_observerBoulder, GpsInView(2:4, i));
    GpsInView(2:4, i) = [az; el; range];
end

azEl_satsT = GpsInView(:, GpsInView(3, :) > 0);
figure
plotAzElCustom(rad2deg(azEl_satsT(2, :)), rad2deg(azEl_satsT(3, :)), azEl_satsT(1, :))
title("GPS satellites over Boulder, CO 1200 26AUG24")

%% 2-2
nineAmMdt = 140418; % week seconds

GpsInView = sp3(sp3(:,8) == 1 & sp3(:, 2) > nineAmMdt & sp3(:, 2) < noonMdt, 2:6)';
GpsInView = [GpsInView(1:2, :); GpsInView(3:5, :)*1000]; % km -> m
for i = 1:size(GpsInView, 2)
    [az, el, range] = computeAzElRange(R_observerBoulder, GpsInView(3:5, i));
    GpsInView(3:5, i) = [az; el; range];
end
GpsInView = GpsInView(:, GpsInView(4, :) > 0);

figure
plotAzElCustom(rad2deg(GpsInView(3, :)), rad2deg(GpsInView(4, :)), GpsInView(2, :))
title("GPS satellites over Boulder, CO 0900-1200 26AUG24")

function plotGPSDifference(DifferenceSeries, title)
    Utilities.multiSeries(DifferenceSeries, ["X" "Y" "Z" "Clock"], ["Time of Ephemeris (s)" "Difference (m)"], title)
end
