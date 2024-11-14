clear; clc; close all;
addpath("..");
GU = GnssUtilities;

%% 1-1
f_sampling = 50e6;               % Sampling frequency: 50 MHz
T_sampling = 1 / f_sampling;              % Sampling period
samplingDuration = 10e-3;        % Duration: 10 ms
t_sampling = 0:T_sampling:samplingDuration;    % Time vector from 0 to 10 ms (exclusive of endpoint)
nfft = length(t_sampling);            % FFT length (same as the number of samples)
df = f_sampling/nfft;            % Frequency resolution
f = df*(0:nfft-1);

%% 2-1
f_signal = 40e3;            % Frequency of sine wave (40 kHz)
A = 1;              % Amplitude (1 V)

signal = A * sin(2 * pi * f_signal * t_sampling);
figure; oscilloscopeViewA1(t_sampling, signal, [0, 1e-3], [0, 1e-4], 'sine 40 kHz');

ps = powerSpectrum(signal);
figure; spectrumAnalyzerView(f, ps, [0, f_sampling/2], [0, 100e3], 'sine 40 kHz');

%% 3-1
% squareSignal = square(2 * pi * f_signal * t);
% figure; oscilloscopeView(t, squareSignal, [0, 1e-3], [0, 1e-4], 'square 40 kHz');

% ps_square = powerSpectrum(squareSignal, nfft);
% figure; spectrumAnalyzerView(f, ps_square, [0, f_sampling/2], [0, 100e3], 'square 40 kHz');

%% 4-1
[~, G1] = GU.generateCA([1, 2], 1023);
t_G1 = 0:length(G1)-1;
t_G1 = t_G1 * 1/1.023e6;

G1Sampled = sampleSignal(G1, t_G1, t_sampling);
G1HwOutput = [t_sampling(485:495)*10^6; G1Sampled(485:495)]

%% 4-2
figure; oscilloscopeViewBinary(t_sampling, G1Sampled, [0, 1e-3], [0, 1e-4], 'G1');
ps_G1 = powerSpectrum(G1Sampled);
figure; spectrumAnalyzerView(f, ps_G1, [0, f_sampling/2], [0, 2024], 'G1');

%% 5-1
phaseSelectors_9 = [3 10];
CA_9 = GU.generateCA(phaseSelectors_9, 1023);
t_CA9 = 0:length(CA_9)-1;
t_CA9 = t_CA9 * 1/1.023e6;
CA_9Sampled = sampleSignal(CA_9, t_CA9, t_sampling);
CA_9HwOutput = [t_sampling(485:495)*10^6; CA_9Sampled(485:495)]

%% 5-2
figure; oscilloscopeViewBinary(t_sampling, CA_9Sampled, [0, 1e-3], [0, 1e-4], 'CA 9');
ps_CA9 = powerSpectrum(CA_9Sampled);
figure; spectrumAnalyzerView(f, ps_CA9, [0, f_sampling/2], [0, 2024], 'CA 9');

%% 6-1
f_carrier = 5.115e6;  % Carrier frequency: 5.115 MHz
carrierSignal = cos(2 * pi * f_carrier * t_sampling);

bpskModulatedSignal = bpskModulate(G1Sampled, carrierSignal);

figure; oscilloscopeViewA1(t_sampling, bpskModulatedSignal, [0, 1e-3], [2.3e-5, 2.4e-5], 'BPSK Modulated Signal');
ps_bpsk = powerSpectrum(bpskModulatedSignal);
figure; spectrumAnalyzerView(f, ps_bpsk, [0, f_sampling/2], [4e6, 6.5e6], 'BPSK Modulated Signal');

%% 6-4
f_carrierDoubled = f_carrier*2;  % Carrier frequency: 10.23 MHz
carrierSignalDoubled = cos(2 * pi * f_carrierDoubled * t_sampling);
bpskModulatedSignalDoubled = bpskModulate(G1Sampled, carrierSignalDoubled);

figure; oscilloscopeViewA1(t_sampling, bpskModulatedSignalDoubled, [0, 1e-3], [2.3e-5, 2.4e-5], 'BPSK Modulated Signal (Carrier Doubled)');
ps_bpskDoubled = powerSpectrum(bpskModulatedSignalDoubled);
figure; spectrumAnalyzerView(f, ps_bpskDoubled, [0, f_sampling/2], [4e6, 6.5e6], 'BPSK Modulated Signal (Carrier Doubled)');

%% 7-1
noise = 1 * randn(size(bpskModulatedSignal));  % Generate white noise with standard deviation 1V
bpskModulatedSignalNoisy = bpskModulatedSignal + noise;  % Add noise to the BPSK modulated signal

figure; oscilloscopeView(t_sampling, bpskModulatedSignalNoisy, [0, 1e-3], [2.3e-5, 2.4e-5], 'BPSK Modulated Signal with Noise');
ps_bpskNoisy = powerSpectrum(bpskModulatedSignalNoisy);
figure; spectrumAnalyzerView(f, ps_bpskNoisy, [0, f_sampling/2], [4e6, 6.5e6], 'BPSK Modulated Signal with Noise');

%% 7-3
bpskUnmodulatedSignal = bpskModulate(G1Sampled, bpskModulatedSignalNoisy);
figure; oscilloscopeView(t_sampling, bpskUnmodulatedSignal, [0, 1e-3], [2.3e-5, 2.4e-5], 'BPSK Unmodulated Signal');
ps_bpskDemodulated = powerSpectrum(bpskUnmodulatedSignal);
figure; spectrumAnalyzerView(f, ps_bpskDemodulated, [0, f_sampling/2], [5e6, 5.2e6], 'BPSK Unmodulated Signal');

function ps = powerSpectrum(signal)
    ps = 20*log10(sqrt(2)*abs(fft(signal)/length(signal)));
end

function oscilloscopeView(t, signal, xlim1, xlim2, title)
    setPosition;

    subplot(1,2,1); plot(t, signal); xlabel('Time (s)'); ylabel('Amplitude (V)'); grid on; xlim(xlim1);
    subplot(1,2,2); plot(t, signal); xlabel('Time (s)'); ylabel('Amplitude (V)'); grid on; xlim(xlim2);
    sgtitle(['Time-Domain Signal (Oscilloscope View) ' title]);
    grid on;
end

function oscilloscopeViewA1(t, signal, xlim1, xlim2, title)
    setPosition;

    subplot(1,2,1); plot(t, signal); xlabel('Time (s)'); ylabel('Amplitude (V)'); grid on; xlim(xlim1); ylim([-1.5 1.5]);
    subplot(1,2,2); plot(t, signal); xlabel('Time (s)'); ylabel('Amplitude (V)'); grid on; xlim(xlim2); ylim([-1.5 1.5]);
    sgtitle(['Time-Domain Signal (Oscilloscope View) ' title]);
    grid on;
end

function oscilloscopeViewBinary(t, signal, xlim1, xlim2, title)
    setPosition;
    binaryAxis = [-.5 1.5];

    subplot(1,2,1); stairs(t, signal); xlabel('Time (s)'); ylabel('Amplitude (V)'); grid on; xlim(xlim1); ylim(binaryAxis);
    subplot(1,2,2); stairs(t, signal); xlabel('Time (s)'); ylabel('Amplitude (V)'); grid on; xlim(xlim2); ylim(binaryAxis);
    sgtitle(['Time-Domain Signal (Oscilloscope View) ' title]);
    grid on;
end

function spectrumAnalyzerView(f, ps, xlim1, xlim2, title)
    setPosition;

    subplot(1,2,1); plot(f, ps, '.-'); xlabel('Frequency (Hz)'); ylabel('Power (dB)'); grid on; xlim(xlim1);
    subplot(1,2,2); plot(f, ps, '.-'); xlabel('Frequency (Hz)'); ylabel('Power (dB)'); grid on; xlim(xlim2);
    sgtitle(['Frequency-Domain Plot (Spectrum Analyzer view) ' title]);
    grid on;
end

function setPosition
    plotPos = [2 2 11 6];
    set(gcf,'unit','inches','position',plotPos);
end

function sampledSignal = sampleSignal(signal, t, t_sample)
    sampledSignal = zeros(1, length(t_sample));
    for i = 1:length(t_sample)
        idx = find(t <= t_sample(i), 1, 'last');
        
        if ~isempty(idx)
            sampledSignal(i) = signal(idx);
        end
    end
end

function bpskModulatedSignal = bpskModulate(binarySignal, carrierSignal)
    signal = 2 * binarySignal - 1;
    bpskModulatedSignal = signal .* carrierSignal;
end