function [aerOpt] = getAerOpt(mTime, height, wavelengths, varargin)
% GETAEROPT get the profiles of aerosol optical properties.
%
% USAGE:
%    [aerOpt] = getAerOpt(height)
%
% INPUTS:
%    mTime: numeric
%        measurement time.
%    height: numeric
%        range bins in vertical. (m)
%    wavelengths: numeric
%        wavelength. (nm)
%
% aerOptModel: char
%    aerosol optical properties calculation model.
%        - clean-urban
%        - moderate-urban
%        - polluted-urban
%        - single-spherical
%        - dust
%        - marine
%        - earlinet
%        - cams
%        - cma: from real-time measurements of cma lidar network.
%        - local-file: from netCDF4 file.
%    aerOptFilepath: datenum
%
% OUTPUTS:
%    aerOpt: struct
%        time: numeric
%        heightArr: numeric
%            vertical range bins. (m)
%        wavelength: numeric
%            wavelength array. (nm)
%        ae: array
%            Angstroem exponent.
%        extinction: matrix (time x height x wavelength)
%            extinction coefficient. (m-1)
%        backscatter: matrix (time x height x wavelength)
%            backscatter coefficient. (m-1 sr-1)
%        depolarization: matrix (time x height x wavelength)
%            particle depolarization ratio.
%
% HISTORY:
%    2024-03-16: first edition by Zhenping
%    2024-08-08: add 'local-file' for input
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'mTime', @isnumeric);
addRequired(p, 'height', @isnumeric);
addRequired(p, 'wavelengths', @isnumeric);
addParameter(p, 'aerOptModel', 'simple', @ischar);
addParameter(p, 'aerOptFilepath', '', @ischar);

parse(p, mTime, height, wavelengths, varargin{:});
%% Aerosol Optical Properties Calculation
switch p.Results.aerOptModel
case 'simple'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = ones(1, length(aerOptRaw.time));

    for iH = 1:length(aerOptRaw.heightArr)
        for iWL = 1:length(aerOptRaw.wavelength)
            if (aerOptRaw.heightArr(iH) <= 2000)
                aerOptRaw.extinction(:, iH, iWL) = (aerOptRaw.wavelength(iWL) / 532) .^ (-aerOptRaw.AE) * 100 * 1e-6;
                aerOptRaw.backscatter(:, iH, iWL) = aerOptRaw.extinction(:, iH, iWL) / 50;
                aerOptRaw.depolarization(:, iH, iWL) = 0.1;
            else
                aerOptRaw.extinction(:, iH, iWL) = 0;
                aerOptRaw.backscatter(:, iH, iWL) = 0;
                aerOptRaw.depolarization(:, iH, iWL) = 0.1;
            end
        end
    end

case 'clean-urban'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = ones;

case 'moderate-urban'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = 1.3 * ones(1, length(aerOptRaw.time));

    for iH = 1:length(aerOptRaw.heightArr)
        for iWL = 1:length(aerOptRaw.wavelength)
            if (aerOptRaw.heightArr(iH) <= 2000)
                aerOptRaw.extinction(:, iH, iWL) = (aerOptRaw.wavelength(iWL) / 532) .^ (-aerOptRaw.AE) * 300 * 1e-6;
                aerOptRaw.backscatter(:, iH, iWL) = aerOptRaw.extinction(:, iH, iWL) / 50;
                aerOptRaw.depolarization(:, iH, iWL) = 0.1;
            else
                aerOptRaw.extinction(:, iH, iWL) = 0;
                aerOptRaw.backscatter(:, iH, iWL) = 0;
                aerOptRaw.depolarization(:, iH, iWL) = 0.1;
            end
        end
    end

case 'polluted-urban'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = 1.3 * ones(1, length(aerOptRaw.time));

    for iH = 1:length(aerOptRaw.heightArr)
        for iWL = 1:length(aerOptRaw.wavelength)
            if (aerOptRaw.heightArr(iH) <= 1500)
                aerOptRaw.extinction(:, iH, iWL) = (aerOptRaw.wavelength(iWL) / 532) .^ (-aerOptRaw.AE) * (2000 / 1.5) * 1e-6;
                aerOptRaw.backscatter(:, iH, iWL) = aerOptRaw.extinction(:, iH, iWL) / 62;
                aerOptRaw.depolarization(:, iH, iWL) = 0.1;
            else
                aerOptRaw.extinction(:, iH, iWL) = 0;
                aerOptRaw.backscatter(:, iH, iWL) = 0;
                aerOptRaw.depolarization(:, iH, iWL) = 0.1;
            end
        end
    end

case 'single-spherical'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = 1.3 * ones(1, length(aerOptRaw.time));

case 'dust'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = 1.3 * ones(1, length(aerOptRaw.time));

case 'marine'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = 1.3 * ones(1, length(aerOptRaw.time));

case 'earlinet'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = 1.3 * ones(1, length(aerOptRaw.time));

case 'cams'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = 1.3 * ones(1, length(aerOptRaw.time));

case 'cma'
    aerOptRaw = struct();
    aerOptRaw.time = mTime;
    aerOptRaw.heightArr = 0:30:40000;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = NaN(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.AE = 1.3 * ones(1, length(aerOptRaw.time));

case 'local-file'
    aerOptIn = struct();
    aerOptIn.time = unix_timestamp_2_datenum(ncread(p.Results.aerOptFilepath, 'time'));
    aerOptIn.heightArr = ncread(p.Results.aerOptFilepath, 'height');
    aerOptIn.wavelength = ncread(p.Results.aerOptFilepath, 'wavelength');
    aerOptIn.extinction = ncread(p.Results.aerOptFilepath, 'extinction');
    aerOptIn.backscatter = ncread(p.Results.aerOptFilepath, 'backscatter');
    aerOptIn.depolarization = ncread(p.Results.aerOptFilepath, 'depolarization_ratio');
    if length(aerOptIn.wavelength) < 2
        error('Insufficient aerosol optical properties. At least two different wavelengths are required.');
    end

    aerOptRaw = struct();
    aerOptRaw.time = aerOptIn.time;
    aerOptRaw.heightArr = aerOptIn.heightArr;
    aerOptRaw.wavelength = wavelengths;
    aerOptRaw.extinction = zeros(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.backscatter = zeros(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));
    aerOptRaw.depolarization = zeros(length(aerOptRaw.time), length(aerOptRaw.heightArr), length(aerOptRaw.wavelength));

    for iWL = 1:length(aerOptRaw.wavelength)
        diffArr = abs(aerOptIn.wavelength - aerOptRaw.wavelength(iWL));
        [~, idx] = sort(diffArr);
        
        extAE = ones(length(aerOptRaw.time), length(aerOptRaw.heightArr));
        bscAE = ones(length(aerOptRaw.time), length(aerOptRaw.heightArr));
        extInt = zeros(length(aerOptRaw.time), length(aerOptRaw.heightArr));
        bscInt = zeros(length(aerOptRaw.time), length(aerOptRaw.heightArr));

        extinction1 = aerOptIn.extinction(:, :, idx(1));
        extinction2 = aerOptIn.extinction(:, :, idx(2));
        backscatter1 = aerOptIn.backscatter(:, :, idx(1));
        backscatter2 = aerOptIn.backscatter(:, :, idx(2));
        isExtValid = (extinction1 > 1e-9) & (extinction2 > 1e-9);
        extAE(isExtValid) = - (log(extinction1(isExtValid)) - log(extinction2(isExtValid))) ./ (log(aerOptIn.wavelength(idx(1))) - log(aerOptIn.wavelength(idx(2))));
        extInt(isExtValid) = extinction1(isExtValid) .* (aerOptRaw.wavelength(iWL) / aerOptIn.wavelength(idx(1))) .^ (-extAE(isExtValid));

        isBscValid = (backscatter1 > 1e-9) & (backscatter2 > 1e-9);
        bscAE(isBscValid) = - (log(backscatter1(isBscValid)) - log(backscatter2(isBscValid))) ./ (log(aerOptIn.wavelength(idx(1))) - log(aerOptIn.wavelength(idx(2))));
        bscInt(isBscValid) = backscatter1(isBscValid) .* (aerOptRaw.wavelength(iWL) / aerOptIn.wavelength(idx(1))) .^ (-bscAE(isBscValid));

        depolInt = aerOptIn.depolarization(:, :, idx(1));

        aerOptRaw.extinction(:, :, iWL) = extInt;
        aerOptRaw.backscatter(:, :, iWL) = bscInt;
        aerOptRaw.depolarization(:, :, iWL) = depolInt;

    end

otherwise
    error('Unknown aerosol optical model: %s', p.Results.aerOptModel);
end

%% Interpolation
aerOpt = struct();
aerOpt.time = aerOptRaw.time;
aerOpt.heightArr = height;
aerOpt.wavelength = aerOptRaw.wavelength;
aerOpt.extinction = NaN(length(aerOpt.time), length(aerOpt.heightArr), length(aerOpt.wavelength));
aerOpt.backscatter = NaN(length(aerOpt.time), length(aerOpt.heightArr), length(aerOpt.wavelength));
aerOpt.depolarization = NaN(length(aerOpt.time), length(aerOpt.heightArr), length(aerOpt.wavelength));

if length(aerOpt.time) > 1
    [TIME, HEIGHT] = meshgrid(aerOpt.heightArr, aerOpt.time);
    [TIMERAW, HEIGHTRAW] = meshgrid(aerOptRaw.heightArr, aerOptRaw.time);
    for iWL = 1:length(aerOpt.wavelength)
        aerOpt.extinction(:, :, iWL) = interp2(double(TIMERAW), double(HEIGHTRAW), double(aerOptRaw.extinction(:, :, iWL)), double(TIME), double(HEIGHT), 'nearest');
        aerOpt.backscatter(:, :, iWL) = interp2(double(TIMERAW), double(HEIGHTRAW), double(aerOptRaw.backscatter(:, :, iWL)),double(TIME), double(HEIGHT), 'nearest');
        aerOpt.depolarization(:, :, iWL) = interp2(double(TIMERAW), double(HEIGHTRAW), double(aerOptRaw.depolarization(:, :, iWL)),double(TIME), double(HEIGHT), 'nearest');
    end
else
    for iWL = 1:length(aerOpt.wavelength)
        aerOpt.extinction(1, :, iWL) = interp1(aerOptRaw.heightArr, squeeze(aerOptRaw.extinction(1, :, iWL)), aerOpt.heightArr, 'nearest');
        aerOpt.backscatter(1, :, iWL) = interp1(aerOptRaw.heightArr, squeeze(aerOptRaw.backscatter(1, :, iWL)), aerOpt.heightArr, 'nearest');
        aerOpt.depolarization(1, :, iWL) = interp1(aerOptRaw.heightArr, squeeze(aerOptRaw.depolarization(1, :, iWL)), aerOpt.heightArr, 'nearest');
    end
end

end