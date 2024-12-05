clear;clc;close all;warning off

% configurations
filename = 'cyg07_raw_if_s20231201_162558_e20231201_162659_data_chn1_first_10sec.dat';
dataFile.sampling = 16.0362e6;
dataFile.dataType = 'bit4';
dataFile.IF = 3.8724e6;

% read data, all
fid = fopen(filename, 'r');
samples_10sec = fread(fid, dataFile.sampling*10, dataFile.dataType)';%10sec

% data visualization of the first 100ms
fs = dataFile.sampling; 
tt = 0:1/fs:1e-1-1/fs; % 100 ms
samples_100ms = samples_10sec(1:length(tt)); clear samples_10sec;
% oscilloscope
figure(201); 
pos = [2 2 11 4];
set(gcf,'unit','inches','position',pos);
subplot(1,2,1);plot(tt,samples_100ms,'.-');xlabel('time (s)');ylabel('amplitude (v)');grid on;xlim()
subplot(1,2,2);plot(tt,samples_100ms,'.-');xlabel('time (s)');ylabel('amplitude (v)');grid on;xlim([0,2e-6])
% spectrum analyzer 
df = fs/length(tt);
f = df*(0:length(tt)-1);
ps = 20*log10(abs(fft(samples_100ms)/length(tt))); 
figure(202); 
pos = [2 2 5.5 4];
set(gcf,'unit','inches','position',pos);
plot(f,ps,'.-');xlabel('freq (Hz)');ylabel('power (dB)');grid on;xlim([0,fs/2])

% % acquisition of GPS L1 C/A signal, 1-ms
% acq_code = nan(1,32);
% acq_Dop = nan(1,32);
% acq_length = 1; % ms
% load('L1CACode.mat')
% 
% acq_2D_corr_record = nan(1,32);
% for prn = 16
%     acq_samples = samples(1:acq_length*1e-3*param.sampling);
%     [acq_2D_corr,search_freqs] = ...
%         acq_corr_fft(acq_samples,param,code(prn,:));
%     acq_2D_corr_record(prn) = max(max(acq_2D_corr));
% end

