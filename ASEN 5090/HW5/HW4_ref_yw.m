clear;clc;close all;warning off;
% reference script for ASEN5090 HW4
% Yang Wang/2024-10
%
addpath UNB3m
c   = 299792458;    % SI speed of light, m/s
omegaE = 7.2921151467e-5;
f1 = 1575.42e6;
f2 = 1227.6e6;

% receiver location
% nist_ecef = [-1288398.575 -4721696.933 4078625.347]';
selected_site_ecef = [-1288398.574, -4721696.936, 4078625.349];

% read ephemeris
navfilename='brdc2390.24n';
[gps_ephem,ionoparams]=read_GPSbroadcast(navfilename);

% sp3
sp3filename='COD0MGXFIN_20242390000_01D_05M_ORB.SP3';
sp3=read_sp3(sp3filename);
sp3(sp3==9.999999999990000e+05)=nan; % set default (useless) values to nan
sp3(:,4:6) = sp3(:,4:6)*1e3;
sp3(:,7) = sp3(:,7)*1e-6*c;

% load rinex file
load('data_rinex_NIST00USA_R_2024_239.mat')

% Q1
figure(101)
dura =  duration(0,0,30);
for selected_prn = [2 16]
    index_prn = data_rinex.GPS.SatelliteID==selected_prn;
    time_rinex = data_rinex.GPS.Time(index_prn);
    C1C = data_rinex.GPS.C1C(index_prn);
    L1C = data_rinex.GPS.L1C(index_prn)*c/f1;
    L2W = data_rinex.GPS.L2W(index_prn)*c/f2;
    sTEC_l = -f1^2*f2^2/40.3/(f1^2-f2^2)*(L2W-L1C);
    iono_L1_l = 40.3*sTEC_l/f1^2;
    % analyze multipath
    % pseudorange carrier phase leveling
    C1C_MP = nan(size(C1C));
    index_jump = unique([1;find(diff(time_rinex)>dura)+1]);
    starts = index_jump(1:end);
    ends = [index_jump(2:end)-1;length(time_rinex)];
    for ii = 1:length(starts)
        C1C_MP(starts(ii):ends(ii)) = C1C(starts(ii):ends(ii)) - L1C(starts(ii):ends(ii)) -...
            2*iono_L1_l(starts(ii):ends(ii));
        C1C_MP_temp = C1C_MP(starts(ii):ends(ii));
        C1C_MP(starts(ii):ends(ii)) = C1C_MP(starts(ii):ends(ii)) ...
            - mean(C1C_MP_temp(~isnan(C1C_MP_temp)));
    end
    plot(time_rinex,C1C_MP,'.'); hold on
    ylim([-3 3]); legend({'PRN2','PRN16'});
end


% Q2
figure(202)
for selected_prn = 8
    index_prn = data_rinex.GPS.SatelliteID==selected_prn;
    time_rinex = data_rinex.GPS.Time(index_prn);
    C1C = data_rinex.GPS.C1C(index_prn);
    C2W = data_rinex.GPS.C2W(index_prn);
    [gps_week_sequence,tow_sequence,gps_week_sequence2]= cal2gps(time_rinex);
    % calculate expected geometry range
    satPos = nan(length(tow_sequence),3);
    satVel = nan(length(tow_sequence),3);
    satClkCorr = nan(size(tow_sequence));
    relCorr = nan(size(tow_sequence));
    tgd = nan(size(tow_sequence));
    Range = nan(size(tow_sequence));
    Az = nan(size(tow_sequence)); El = nan(size(tow_sequence));
    for ij = 1:length(tow_sequence)
        [health0,satPos0,satVel0,satClkCorr0,relCorr0,tgd0] = ...
            eph2pvt(gps_ephem,[gps_week_sequence(ij),tow_sequence(ij)],selected_prn);
        [Az(ij),El(ij),Range(ij)] = compute_azelrange(selected_site_ecef, satPos0);
        range = norm(selected_site_ecef' - satPos0');
        for ik = 1:3
            tt = tow_sequence(ij)-range/c;
            [health0,satPos0,satVel0,satClkCorr0,relCorr0,tgd0] = ...
                eph2pvt(gps_ephem,[gps_week_sequence(ij),tt],selected_prn);
            phi = (tow_sequence(ij)-tt)*omegaE;
            satPos0 = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1; ]*satPos0';
            range = norm(selected_site_ecef' - satPos0);
        end
        satPos(ij,:) = satPos0';
        satVel(ij,:) = satVel0';
        satClkCorr(ij) = satClkCorr0;
        relCorr(ij) = relCorr0;
        tgd(ij) = tgd0;
        [Az(ij),El(ij),Range(ij)] = compute_azelrange(selected_site_ecef, satPos0');
    end
    % calculate tropospehric delay
    % unb3m
    doy = day(time_rinex,'dayofyear');
    timeInHours = mod(tow_sequence,24*3600)/60/60; % hr
    nist_lla = ecef2lla(selected_site_ecef);
    RTROP = nan(size(tow_sequence));
    orth_height = nist_lla(3)-geoidheight(nist_lla(1),nist_lla(2));
    for ij = 1:length(tow_sequence)
        [RTROP(ij) HZD(ij) HMF(ij) WZD(ij) WMF(ij)]...
            =UNB3M(nist_lla(1)/180*pi,orth_height,doy(ij)+timeInHours(ij)/24,El(ij)/180*pi);
    end
    % calculate ionospheric delay
    sTEC_c = f1^2*f2^2/40.3/(f1^2-f2^2)*(C2W-C1C);
    iono_L1_c = 40.3*sTEC_c/f1^2;
    C_if = C1C-iono_L1_c;
    % residual 
    rho_res = C_if-Range-satClkCorr-relCorr-RTROP;
    plot(time_rinex,rho_res,'.');hold on
%     ylim([50 80])
end

if selected_prn==8
    figure(201)
    subplot(3,2,1); plot(time_rinex,Range,'.'); ylabel('Expected range (m)');
    subplot(3,2,2); plot(time_rinex,satClkCorr,'.'); ylabel('Sat Clk Bias (m)');
    subplot(3,2,3); plot(time_rinex,relCorr,'.'); ylabel('Relat. Corr (m)');
    subplot(3,2,4); plot(time_rinex,iono_L1_c,'.'); ylabel('Iono. delay (m)');
    subplot(3,2,5); plot(time_rinex,RTROP,'.'); ylabel('Tropo. delay (m)'); ylim([0 50]);
    format long
    disp([Range(1) Range(end);satClkCorr(1) satClkCorr(end);relCorr(1) relCorr(end);...
        iono_L1_c(1) iono_L1_c(end);RTROP(1) RTROP(end);])
end


% Q4=3
figure(301)
for selected_prn = 8
    index_prn = data_rinex.GPS.SatelliteID==selected_prn;
    time_rinex = data_rinex.GPS.Time(index_prn);
    D1C = data_rinex.GPS.D1C(index_prn);
    [gps_week_sequence,tow_sequence,gps_week_sequence2]= cal2gps(time_rinex);
    % calculate expected geometry range
    for ij = 1:length(tow_sequence)
        [health0,satPos0,satVel0,satClkCorr0,relCorr0,tgd0] = ...
            eph2pvt(gps_ephem,[gps_week_sequence(ij),tow_sequence(ij)],selected_prn);
        range0 = norm(selected_site_ecef' - satPos0');
        for ik = 1:3
            tt = tow_sequence(ij)-range/c;
            [health0,satPos1,satVel1,satClkCorr0,relCorr0,tgd0] = ...
                eph2pvt(gps_ephem,[gps_week_sequence(ij),tt],selected_prn);
            phi = (tow_sequence(ij)-tt)*omegaE;
            satPos1 = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1; ]*satPos1';
            satVel1 = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1; ]*satVel1';
            range = norm(selected_site_ecef' - satPos1);
        end
        los_GPS_nist = (selected_site_ecef'-satPos1)/range;
        Doppler0(ij) = satVel0*los_GPS_nist/c*f1;
        Doppler(ij) = satVel1'*los_GPS_nist/c*f1;
    end
    plot(time_rinex,(Doppler-D1C')/f1*c,'.'); ylabel('Doppler (m/s)'); hold on;
    plot(time_rinex,(Doppler0-D1C')/f1*c,'.'); ylabel('Doppler (m/s)')
    plot(time_rinex,(Doppler0-Doppler)/f1*c,'.'); ylabel('Doppler (m/s)')
% 
% disp(Doppler(ij))
% disp(Doppler(ij)/f1*c)
% disp(satVel1)
% disp(los_GPS_nist)

end
    







