function [fileName] = write_CMA_L0(savePath, stationNumber, producerID, deviceNumber, longitude, latitude, elevation, startTime, endTime, elevationAngle, wavelength1, wavelength2, wavelength3, numChannel, digitizer, recWavelength, recType, resolution, overlap, nBins, lidarSignal, varargin)
% WRITE_CMA_L0 write multi-param lidar data based on CMA BIN format.
%
% USAGE:
%    [fileName] = write_CMA_L0(savePath, stationNumber, producerID, deviceNumber, longitude, latitude, elevation, startTime, endTime, elevationAngle, wavelength1, wavelength2, wavelength3, numChannel, digitizer)
%
% INPUTS:
%    savePath: char
%        save path.
%    stationNumber: char
%        five digit station number (e.g, '57511')
%    producerID: char
%        producer ID (e.g., 'ILYJ1')
%    deviceNumber: numeric
%        device number.
%    longitude: numeric
%        longitude (degree)
%    latitude: numeric
%        latitude (degree)
%    elevation: numeric
%        elevation (m)
%    startTime: numeric
%        starttime of the measurement (matlab datenum)
%    endTime: numeric
%        stop time of the measurement (matlab datenum)
%    elevationAngle: numeric
%        elevation angle (degree)
%    wavelength1: numeric
%        emission wavelength 1 (nm)
%    wavelength2: numeric
%        emission wavelength 2 (nm)
%    wavelength3: numeric
%        emission wavelength 3 (nm)
%    numChannel: numeric
%        number of receiving channels.
%    digitizer: array
%        digitizer mode for each channel. (0: AD; 1: PC; 2: Combined)
%    recWavelength: array
%        central wavelengths for each channel. (nm)
%    recType: array
%        receiving type for each channel. (0: elastic; 1: P; 2: S; 3: Raman)
%    resolution: array
%        resolution of one range bin for each channel. (m)
%    overlap: array
%        height with complete overlap. (m)
%    nBins: array
%        range bins for each channel.
%    lidarSignal: matrix (channel x bin)
%
% KEYWORDS:
%    debug: logical
%
% OUTPUTS:
%    fileName: char
%        filename of the output file.
%
% REFERENCES:
%    ./docs/附录A 温湿度气溶胶激光雷达数据格式规范.docx
%
% HISTORY:
%    2024-06-10: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'savePath', @ischar);
addRequired(p, 'stationNumber', @ischar);
addRequired(p, 'producerID', @ischar);
addRequired(p, 'deviceNumber', @isnumeric);
addRequired(p, 'longitude', @isnumeric);
addRequired(p, 'latitude', @isnumeric);
addRequired(p, 'elevation', @isnumeric);
addRequired(p, 'startTime', @isnumeric);
addRequired(p, 'endTime', @isnumeric);
addRequired(p, 'elevationAngle', @isnumeric);
addRequired(p, 'wavelength1', @isnumeric);
addRequired(p, 'wavelength2', @isnumeric);
addRequired(p, 'wavelength3', @isnumeric);
addRequired(p, 'numChannel', @isnumeric);
addRequired(p, 'digitizer', @isnumeric);
addRequired(p, 'recWavelength', @isnumeric);
addRequired(p, 'recType', @isnumeric);
addRequired(p, 'resolution', @isnumeric);
addRequired(p, 'overlap', @isnumeric);
addRequired(p, 'nBins', @isnumeric);
addRequired(p, 'lidarSignal', @isnumeric);
addParameter(p, 'debug', false, @islogical);

parse(p, savePath, stationNumber, producerID, deviceNumber, longitude, latitude, elevation, startTime, endTime, elevationAngle, wavelength1, wavelength2, wavelength3, numChannel, digitizer, recWavelength, recType, resolution, overlap, nBins, lidarSignal, varargin{:});

%% Write data
fileName = sprintf('Z_RADA_I_%s_%s_O_MLIDAR_%s_L0.BIN', stationNumber, datestr(startTime, 'yyyymmddHHMMSS'), producerID);
fid = fopen(fullfile(savePath, fileName), 'w');

fwrite(fid, zeros(1, 7), 'short');   % reserve
fwrite(fid, 0, 'ushort');   % data header
fwrite(fid, 1, 'ushort');   % version
fwrite(fid, deviceNumber, 'uint');   % device number
fwrite(fid, longitude * 1e4, 'uint');   % longitude
fwrite(fid, latitude * 1e4, 'uint');   % latitude
fwrite(fid, elevation * 1e2, 'uint');   % elevation
fwrite(fid, 0, 'ushort');   % reserve
fwrite(fid, 1, 'ushort');   % detection mode

startSeconds = mod(startTime, 1) * 3600 * 24;   % start time in seconds
endSeconds = mod(endTime, 1) * 3600 * 24;   % end time in seconds
startDate = floor(startTime - datenum(1970, 1, 1));   % julian date
fwrite(fid, startSeconds, 'uint');
fwrite(fid, endSeconds, 'uint');
fwrite(fid, startDate, 'ushort');
fwrite(fid, elevationAngle / (180 / 4096) * 8, 'ushort');   % elevation angle
fwrite(fid, 0, 'ushort');   % reserve
fwrite(fid, wavelength1, 'ushort');   % wavelength 1
fwrite(fid, wavelength2, 'ushort');   % wavelength 2
fwrite(fid, wavelength3, 'ushort');   % wavelength 3
fwrite(fid, numChannel, 'ushort');   % number of channels

% channel header
for iCh = 1:numChannel
    fwrite(fid, iCh, 'ushort');   % channel index
    fwrite(fid, bitor(bitshift(digitizer(iCh), 14), recWavelength(iCh)), 'ushort');
    fwrite(fid, recType(iCh), 'ushort');   % receiving  channel type (0: elastic; 1: P; 2: S; 3: Raman)
    fwrite(fid, resolution(iCh) * 100, 'ushort');   % resolution
    fwrite(fid, overlap(iCh) * 10, 'ushort');   % overlap height
    fwrite(fid, 31 + 17 * numChannel + 1 + (iCh - 1) * nBins(1), 'uint');
    fwrite(fid, nBins(iCh), 'ushort');
end

% write data
for iCh = 1:numChannel
    fwrite(fid, lidarSignal(iCh, :), 'float');
end

fclose(fid);

end