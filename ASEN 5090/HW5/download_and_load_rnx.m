% download and read hourly rinex 3 files from cddis for a day
% https://cddis.nasa.gov/archive/gnss/data/hourly/
% Y Wang (2023/09/27)
%
clear;clc;close all;warning off;

% specify doy and hh in the for loops
% also specify the IGS site name in the script below
% e.g., here is NIST00USA_R_...
% edit the directory for CRX2RNX.exe 

site = 'NIST00USA_R_';
year = '2024';
mkdir('data')
data_rinex = [];

for doy = 239
    for ij = 0:23

        % downloading cmd line
        filename = ['NIST00USA_R_' year num2str(doy,'%03d') num2str(ij,'%02d') ...
            '00_01H_30S_MO.crx.gz'];
        cmd_line = ['curl -u anonymous:test@gmail.com -O --ftp-ssl ' ...
            '--max-time 5 ' ...
            'ftp://gdc.cddis.eosdis.nasa.gov/gnss/data/hourly/' year '/' ...
            num2str(doy,'%03d') '/' num2str(ij,'%02d') '/' ...
            filename];
        system(cmd_line)

        % lots of sites use compressed rinex format (crx)
        % this needs to be unzipped with 'crx2rnx.exe'
        % https://terras.gsi.go.jp/ja/crx2rnx.html
        gunzip(filename)
        cmd_line = ['CRX2RNX.exe ' filename(1:end-3)];
        system(cmd_line)

        % read hourly data into a day
        data_temp = rinexread([filename(1:end-6) 'rnx']);
        if isempty(data_rinex)
            varaible_names = fieldnames(data_temp);
            data_rinex = data_temp;
        else
            for ii = 1:size(varaible_names)
                data_rinex.(string(varaible_names(ii))) = ...
                    [data_rinex.(string(varaible_names(ii))); ...
                    data_temp.(string(varaible_names(ii)))];
            end
        end

        % move files to 'data'
        movefile(filename, 'data')
        movefile(filename(1:end-3), 'data')
        movefile([filename(1:end-6) 'rnx'], 'data')
    end
end

clear data_temp varaible_names
save(['data_rinex_' site year '_' num2str(doy,'%03d')], 'data_rinex')

