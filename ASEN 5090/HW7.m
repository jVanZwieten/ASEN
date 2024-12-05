close all; clear; clc; addpath("..")
GU = GnssUtilities;

%% 2
R_earth = CelestialParameters.radius_earth/1000; % km
mu_earth = CelestialParameters.gravityParameter_earth;

alt_Cygnss = 520; % km
r_Cygnss = R_earth + alt_Cygnss
v_Cygnss = sqrt(mu_earth / r_Cygnss)

alt_Gps = 20200; % km
r_Gps = R_earth + alt_Gps
v_Gps = sqrt(mu_earth / r_Gps)

v_GpsLos = v_Gps*r_Cygnss/r_Gps
v_relMax = v_GpsLos + v_Cygnss

c = Utilities.speedOfLight/1000; % km/s
f_L1 = 1575.42e6; % Hz
f_doppler = f_L1*v_relMax/c

% %% 3a
chippingRate_gps = GnssUtilities.gpsChippingRate; % Hz
sampleRate = 16.0362e6; % Hz
samplePeriod = 1/sampleRate; % s
t_sample = 0:samplePeriod:1e-3-samplePeriod; % 1 ms

filename = 'cyg07_raw_if_s20231201_162558_e20231201_162659_data_chn1_first_10sec.dat';
dataType = 'bit4';
fid = fopen(filename, 'r');
CygnssSamples_10sec = fread(fid, sampleRate*10, dataType)';

% nfft = length(t_sample);         % FFT length
% df = sampleRate/nfft;            % Frequency resolution
% fDomain = df*(0:nfft-1);

% %% 3b
[~, CA_prn32] = get_CA_code(32);
t_CA = (0:length(CA_prn32)-1)/chippingRate_gps;
CA_prn32Sampled = GU.sampleSignal(CA_prn32, t_CA, t_sample);

% %% 3c
f_intermediate = 3.8724e6; % Hz
f_doppler = 10.5e3; % Hz
f_carrier = f_intermediate + f_doppler;

% %% 3d
referenceSignal = exp(-2*pi*f_carrier*t_sample*j).*CA_prn32Sampled;
% powerSpectrum_32I = GU.powerSpectrum(real(referenceSignal));
% figure; GU.spectrumAnalyzerView(fDomain, powerSpectrum_32I, [0, sampleRate/2]); title('Power Spectrum of Replicated PRN 32 (I)');

% %% 3e
CygnssSamples_1ms = CygnssSamples_10sec(1:length(t_sample));

% tic
[correlation, lags] = cyc_corr_basic(CygnssSamples_1ms, referenceSignal);
% t_1msCorrelation = toc

correlation = abs(correlation);
[maxCorrelation, peakIndex] = max(correlation);
% figure; plot(lags*samplePeriod, correlation); xlabel('t (s)'); ylabel('Cross-Correlation'); title('Correlation, 1ms');

% peakIndexChip = peakIndex*samplePeriod*chippingRate_gps;
% peakChip = mod(peakIndexChip, 1023)
% code_phase = 1023 - peakChip

% %% 3f
% t_sample_5ms = linspace(0, 5e-3, 80180);
% CygnssSamples_5ms = CygnssSamples_10sec(1:length(t_sample_5ms));
% CA_prn32Sampled_5ms = [CA_prn32Sampled CA_prn32Sampled CA_prn32Sampled CA_prn32Sampled CA_prn32Sampled];
% referenceSignal_5ms = exp(-2*pi*f_carrier*t_sample_5ms*j).*CA_prn32Sampled_5ms;

% tic
% [correlation_5ms, lags_5ms] = cyc_corr_basic(CygnssSamples_5ms, referenceSignal_5ms);
% t_5msCorrelation = toc

% correlation_5ms = abs(correlation_5ms);
% [maxCorrelation_5ms, peakIndex_5ms] = max(correlation);

% figure
% plot(t_sample_5ms, correlation_5ms(1:length(t_sample_5ms)));
% xlabel('t (s)');
% ylabel('Cross-Correlation');
% title('Correlation, 5ms');

% figure;
% plot(lags * samplePeriod * chippingRate_gps, correlation);
% xlabel('Chip');
% ylabel('Cross-Correlation');
% title('Correlation 0-1023 chips, 1ms');
% xlim([0 1023]);

% figure;
% plot(lags_5ms * samplePeriod * chippingRate_gps, correlation_5ms);
% xlabel('Chip');
% ylabel('Cross-Correlation');
% title('Correlation 0-1023 chips, 5ms');
% xlim([0 1023]);

% %% 4
% tic
CygnssFFT_1ms = fft(CygnssSamples_1ms);
referenceFFT_1ms = conj(fft(referenceSignal));
correlationFFT_1ms = ifft(CygnssFFT_1ms.*referenceFFT_1ms);
% t_1msCorrelationFFT = toc
correlationFFT_1ms = abs(correlationFFT_1ms);
[maxCorrelationFFT_1ms, peakIndexFFT_1ms] = max(correlationFFT_1ms);

peakIndexChip = peakIndex*samplePeriod*chippingRate_gps;
peakChip = mod(peakIndexChip, 1023);
codePhase_fft = 1023 - peakChip

% figure;
% plot(t_sample, correlationFFT_1ms(1:length(t_sample)));
% xlabel('t (s)');
% ylabel('Cross-Correlation');
% title('Correlation via FFT, 1ms');

%% 5a
deltaF_doppler1ms = 1/(2*1e-3)
deltaF_doppler5ms = 1/(2*5e-3)

%% 5b
[~, CA_prn3] = get_CA_code(3);
t_CA3 = (0:length(CA_prn3)-1)/chippingRate_gps;
CA_prn3Sampled = GU.sampleSignal(CA_prn3, t_CA, t_sample);

[correlationScan_prn3, f_dopplers, t_correlation] = correlationDopplerScan([-45.3e3, 45.3e3], deltaF_doppler1ms, f_intermediate, CA_prn3Sampled, CygnssFFT_1ms, t_sample);

figure;
surf(t_correlation, f_dopplers, correlationScan_prn3);
xlabel('lag (chip)');
ylabel('Doppler Frequency (Hz)');
zlabel('Correlation');
title('Correlation Scan');

%% 5c
correlation_values = correlationScan_prn3(:);
figure;
histogram(correlation_values, 50);
xlabel('Correlation Magnitude');
ylabel('Frequency');
title('Histogram of Correlation Results');

%% 5d
prnsPresentTab = table('Size', [0 3], 'VariableTypes', {'double', 'double', 'double'}, 'VariableNames', {'prn', 'doppler', 'lag'});
for prn = 1:31
    [~, CA_prn] = get_CA_code(prn);
    CA_prnSampled = GU.sampleSignal(CA_prn, t_CA, t_sample);

    [f_doppler, lag] = acquireSignal(prn, [-45.3e3, 45.3e3], deltaF_doppler5ms, f_intermediate, CygnssFFT_1ms, t_sample);

    if ~isnan(f_doppler)
        prnsPresentTab = [prnsPresentTab; {prn, f_doppler, lag}];
    end
end

prnsPresentTab = [prnsPresentTab; {32, 10.5e3, codePhase_fft}]

function [unique_t_chipDomain, correlation_chipDomain] = correlationTime2ChipDomain(correlation, t, chippingRate)
    t_chipDomain = t*chippingRate;
    t_chipDomain = mod(t_chipDomain, 1023);

    [t_chipDomainSorted, sortIdx] = sort(t_chipDomain);
    correlationSorted = correlation(sortIdx);

    [unique_t_chipDomain, ~, idx] = unique(t_chipDomainSorted);
    correlation_chipDomain = accumarray(idx, correlationSorted, [], @mean);
end

function [correlationScan, f_dopplers, t_correlation] = correlationDopplerScan(dopplerRange, deltaF, f_intermediate, CA_sampled, signal_acquisitionFft, t_sample)
    chippingRate = GnssUtilities.gpsChippingRate;

    f_dopplers = dopplerRange(1):deltaF:dopplerRange(2);
    for i = 1:length(f_dopplers)
        f_doppler = f_dopplers(i);
        f_carrier = f_intermediate + f_doppler;
        referenceSignal = exp(-2*pi*f_carrier*t_sample*j).*CA_sampled;
        referenceFFT = conj(fft(referenceSignal));

        correlation_fDopler = abs(ifft(signal_acquisitionFft.*referenceFFT));
        [t_correlation, correlation_fDopler] = correlationTime2ChipDomain(correlation_fDopler, t_sample, chippingRate);

        if ~exist('correlationScan', 'var')
            correlationScan = zeros(length(f_dopplers), length(t_correlation));
        end

        correlationScan(i, :) = correlation_fDopler;
    end
end

function [f_doppler, lag] = acquireSignal(prn, dopplerRange, deltaF, f_intermediate, signal_acquisitionFft, t_sample)
    chippingRate_gps = GnssUtilities.gpsChippingRate;

    [~, CA_prn] = get_CA_code(prn);
    t_CA = (0:length(CA_prn)-1)/chippingRate_gps;
    CA_prnSampled = GnssUtilities.sampleSignal(CA_prn, t_CA, t_sample);

    f_dopplers = dopplerRange(1):deltaF:dopplerRange(2);
    for f_doppler = f_dopplers
        f_carrier = f_intermediate + f_doppler;
        referenceSignal = exp(-2*pi*f_carrier*t_sample*j).*CA_prnSampled;
        referenceFFT = conj(fft(referenceSignal));

        correlation_fDopler = abs(ifft(signal_acquisitionFft.*referenceFFT));
        if any(correlation_fDopler > 2000)
            [~, lagIdx] = max(correlation_fDopler);
            lag = t_sample(lagIdx);
            lag_chip = lag*chippingRate_gps;
            lag_chip = mod(lag_chip, 1023);
            lag = 1023 - lag_chip;
            return;
        end
    end

    f_doppler = NaN;
    lag = NaN;
end