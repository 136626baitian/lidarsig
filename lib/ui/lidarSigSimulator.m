function [lidarSig] = lidarSigSimulator(mTime, lidarCfgFile, varargin)
% LIDARSIGSIMULATOR lidar signal simulation.
%
% USAGE:
%    [lidarSig] = lidarSigSimulator(lidarCfgFile, obsCfgFile)
%
% INPUTS:
%    lidarCfgFile: char
%        lidar configuration file.
%
% KEYWORDS:
%    obsCfgFile: char
%        observation configuration file.
%    savePath: char
%        save path.
%    flagRemoveSignalFiles
%    visible: char
%        'on' or 'off'
%    lidarCfg: struct
%        if not empty, lidarCfg in the signal simulation will be derived from this keyword.
%    obsCfg: struct
%
% OUTPUTS:
%    lidarSig: struct
%       level0.sig
%       level1.altitude
%       level1.range
%       level1.sig
%       level1.bg
%
% EXAMPLE:
%
% HISTORY:
%    2024-04-03: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'mTime', @isnumeric);
addRequired(p, 'lidarCfgFile', @ischar);
addParameter(p, 'obsCfgFile', '', @ischar);
addParameter(p, 'savePath', '', @ischar);
addParameter(p, 'visible', 'off', @ischar);
addParameter(p, 'lidarCfg', struct(), @isstruct);
addParameter(p, 'obsCfg', struct(), @isstruct);

parse(p, mTime, lidarCfgFile, varargin{:});

global LISAR_ENVS;

%% Parameter Definition
lidarCfgDefaultFile = fullfile(LISAR_ENVS.RootDir, 'lib', 'config', 'lidarCfg_default.yml');
obsCfgDefaultFile = fullfile(LISAR_ENVS.RootDir, 'lib', 'config', 'obsCfg_default.yml');
if isempty(p.Results.obsCfgFile)
    obsCfgFile = obsCfgDefaultFile;
else
    obsCfgFile = p.Results.obsCfgFile;
end
const = loadConstants();

%% Read Configuration File
if isempty(fieldnames(p.Results.lidarCfg))
    lidarCfg = load_settings(lidarCfgFile, lidarCfgDefaultFile);
else
    lidarCfg = p.Results.lidarCfg;
end
if isempty(fieldnames(p.Results.obsCfg))
    obsCfg = load_settings(obsCfgFile, obsCfgDefaultFile);
else
    obsCfg = p.Results.obsCfg;
end

nBins = (lidarCfg.generalCfg.binNum - lidarCfg.generalCfg.nPretrigger);
hRes = lidarCfg.generalCfg.binWidth * const.c * 5e-10;
distArr = ((1:nBins) - 0.5) * hRes;
heightArr = distArr * cos(lidarCfg.generalCfg.zenithAngle / 180 * pi);
altArr = heightArr + obsCfg.geoCfg.altitude;

%% Input Checker
nPrf = length(mTime);
deltaTDef = inf;
if nPrf > 1
    deltaTDef = (mTime(2) - mTime(1)) * 24 * 3600;   % temporal interval between adjacent time. (s)
end

accTLidar = lidarCfg.generalCfg.accShots * (1 / lidarCfg.generalCfg.laserRepRate);
if accTLidar > deltaTDef
    error('mTime interval is not enough for accumulating sufficient laser shots. (min: %f s)', accTLidar);
end

%% Read Meteorological Data
meteorData = readMeteorData(mTime, obsCfg.meteorCfg);

% data clean
NPRF = length(meteorData.time);
NMAXLEVELS = 1000;
MAXVALIDPOINTS = inf;
meteorDataClean = struct();
meteorDataClean.time = NaN(NPRF, NMAXLEVELS);
meteorDataClean.altitude = NaN(NPRF, NMAXLEVELS);
meteorDataClean.temperature = NaN(NPRF, NMAXLEVELS);
meteorDataClean.pressure = NaN(NPRF, NMAXLEVELS);
meteorDataClean.water_vapor_mixing_ratio = NaN(NPRF, NMAXLEVELS);
for iPrf = 1:NPRF
    temperature = meteorData.temperature(:, iPrf);
    pressure = meteorData.pressure(:, iPrf);
    water_vapor_mixing_ratio = meteorData.water_vapor_mixing_ratio(:, iPrf);
    altitude = meteorData.altitude(:, iPrf);

    isValidPts = (~ isnan(temperature)) & (~ isnan(pressure));
    meteorDataClean.time(iPrf, 1:sum(isValidPts)) = meteorData.time(iPrf);
    meteorDataClean.altitude(iPrf, 1:sum(isValidPts)) = altitude(isValidPts);
    meteorDataClean.temperature(iPrf, 1:sum(isValidPts)) = temperature(isValidPts);
    meteorDataClean.pressure(iPrf, 1:sum(isValidPts)) = pressure(isValidPts);
    meteorDataClean.water_vapor_mixing_ratio(iPrf, 1:sum(isValidPts)) = water_vapor_mixing_ratio(isValidPts);
    meteorDataClean.water_vapor_mixing_ratio(iPrf, isnan(meteorDataClean.water_vapor_mixing_ratio(iPrf, :))) = 0;

    nLevels = sum(isValidPts);
    if nLevels < MAXVALIDPOINTS
        MAXVALIDPOINTS = nLevels;
    end
end
%这里代码改变
meteorDataClean.time = meteorDataClean.time(:, 1:MAXVALIDPOINTS);
meteorDataClean.altitude = meteorDataClean.altitude(:, 1:MAXVALIDPOINTS);
meteorDataClean.temperature = meteorDataClean.temperature(:, 1:MAXVALIDPOINTS);
meteorDataClean.pressure = meteorDataClean.pressure(:, 1:MAXVALIDPOINTS);
meteorDataClean.water_vapor_mixing_ratio = meteorDataClean.water_vapor_mixing_ratio(:, 1:MAXVALIDPOINTS);

%% Calculate Molecular Scattering
nRotNum = 20;
mBsc = NaN(NPRF, MAXVALIDPOINTS, (10 * nRotNum + 770 + 1) * length(lidarCfg.generalCfg.laserWavelengths));
mBscWL = NaN(NPRF, MAXVALIDPOINTS, (10 * nRotNum + 770 + 1) * length(lidarCfg.generalCfg.laserWavelengths));
mExtEmit = NaN(NPRF, MAXVALIDPOINTS, length(lidarCfg.generalCfg.laserWavelengths));
mExtRecv = NaN(NPRF, MAXVALIDPOINTS, length(lidarCfg.generalCfg.channelLabels));

for iPrf = 1:NPRF
    meteorTEMP = meteorDataClean.temperature(iPrf, :);
    meteorPRES = meteorDataClean.pressure(iPrf, :);
    meteorWVMR = meteorDataClean.water_vapor_mixing_ratio(iPrf, :);

    for iL = 1:MAXVALIDPOINTS
        counter = 1;

        for iWL = 1:length(lidarCfg.generalCfg.laserWavelengths)

            nAirMol = number_density_at_pt(meteorPRES(iL), meteorTEMP(iL) - const.T0, 40, 1);
            nN2 = nAirMol * const.N2_parameters.relative_concentration;
            nO2 = nAirMol * const.O2_parameters.relative_concentration;
            nH2O = meteorWVMR(iL) * 28.9660 * 1e-3 * nAirMol / 18.0160;

            [thisScaSigma, thisScaWL] = getDiMolRamanLines(...
                lidarCfg.generalCfg.laserWavelengths{iWL}, ...
                meteorTEMP(iL) - const.T0, ...
                'showFigure', false, ...
                'molecular', 'N2', ...
                'branch', 'RR', ...
                'rotNumArray', 0:(nRotNum - 1));
            mBsc(iPrf, iL, counter:(counter + length(thisScaSigma) - 1)) = thisScaSigma * nN2;
            mBscWL(iPrf, iL, counter:(counter + length(thisScaSigma) - 1)) = thisScaWL;
            counter = counter + length(thisScaSigma);

            [thisScaSigma, thisScaWL] = getDiMolRamanLines(...
                lidarCfg.generalCfg.laserWavelengths{iWL}, ...
                meteorTEMP(iL) - const.T0, ...
                'showFigure', false, ...
                'molecular', 'O2', ...
                'branch', 'RR', ...
                'rotNumArray', 0:(nRotNum - 1));
            mBsc(iPrf, iL, counter:(counter + length(thisScaSigma) - 1)) = thisScaSigma * nO2;
            mBscWL(iPrf, iL, counter:(counter + length(thisScaSigma) - 1)) = thisScaWL;
            counter = counter + length(thisScaSigma);

            [thisScaSigma, thisScaWL] = getDiMolRamanLines(...
                lidarCfg.generalCfg.laserWavelengths{iWL}, ...
                meteorTEMP(iL) - const.T0, ...
                'showFigure', false, ...
                'molecular', 'N2', ...
                'branch', 'VRR', ...
                'rotNumArray', 0:(nRotNum - 1));
            mBsc(iPrf, iL, counter:(counter + length(thisScaSigma) - 1)) = thisScaSigma * nO2;
            mBscWL(iPrf, iL, counter:(counter + length(thisScaSigma) - 1)) = thisScaWL;
            counter = counter + length(thisScaSigma);

            [thisScaSigma, thisScaWL] = getDiMolRamanLines(...
                lidarCfg.generalCfg.laserWavelengths{iWL}, ...
                meteorTEMP(iL) - const.T0, ...
                'showFigure', false, ...
                'molecular', 'O2', ...
                'branch', 'VRR', ...
                'rotNumArray', 0:(nRotNum - 1));
            mBsc(iPrf, iL, counter:(counter + length(thisScaSigma) - 1)) = thisScaSigma * nO2;
            mBscWL(iPrf, iL, counter:(counter + length(thisScaSigma) - 1)) = thisScaWL;
            counter = counter + length(thisScaSigma);

            mBsc(iPrf, iL, counter) = beta_pi_rayleigh(...
                lidarCfg.generalCfg.laserWavelengths{iWL}, ...
                meteorPRES(iL), meteorTEMP(iL) - const.T0, 380, 40);
            mBscWL(iPrf, iL, counter) = lidarCfg.generalCfg.laserWavelengths{iWL};
            counter = counter + 1;

            [thisWVLines, thisOutWL] = getWVRamanLines(...
                lidarCfg.generalCfg.laserWavelengths{iWL} * 1e-9, meteorTEMP(iL) - const.T0);
            mBsc(iPrf, iL, counter:(counter + length(thisWVLines) - 1)) = thisWVLines * nH2O;
            mBscWL(iPrf, iL, counter:(counter + length(thisWVLines) - 1)) = thisOutWL * 1e9;
            counter = counter + length(thisOutWL);

            mExtEmit(iPrf, iL, iWL) = alpha_rayleigh(...
                lidarCfg.generalCfg.laserWavelengths{iWL}, ...
                meteorPRES(iL), meteorTEMP(iL) - const.T0, 380, 40);
        end

        for iWL = 1:length(lidarCfg.generalCfg.channelLabels)
            chLabel = lidarCfg.generalCfg.channelLabels{iWL};
            mExtRecv(iPrf, iL, iWL) = alpha_rayleigh(...
                lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, ...
                meteorPRES(iL), meteorTEMP(iL) - const.T0, 380, 40);
        end
    end
end

%% molecular properties interpolation 
mBscInterp = NaN(length(mTime), length(heightArr), length(lidarCfg.generalCfg.channelLabels));
mExtEmitInterp = NaN(length(mTime), length(heightArr), length(lidarCfg.generalCfg.laserWavelengths));
mExtRecvInterp = NaN(length(mTime), length(heightArr), length(lidarCfg.generalCfg.channelLabels));
for iPrf = 1:length(mTime)
    % 获取当前时间点的气象数据
    currentAlt = meteorData.altitude(:, iPrf);
    currentTemp = meteorData.temperature(:, iPrf);
    currentPres = meteorData.pressure(:, iPrf);
    currentWVMR = meteorData.water_vapor_mixing_ratio(:, iPrf);

    % 去除无效值
    validIdx = ~isnan(currentTemp) & ~isnan(currentPres);
    currentAlt = currentAlt(validIdx);
    currentTemp = currentTemp(validIdx);
    currentPres = currentPres(validIdx);
    currentWVMR = currentWVMR(validIdx);

    % 对每个高度参数进行一维插值
    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        chLabel = lidarCfg.generalCfg.channelLabels{iCh};
        % 计算分子散射并插值到altArr
        mBscInterp(iPrf, :, iCh) = interp1(currentAlt, mBsc(iPrf, 1:length(currentAlt), iCh), altArr, 'linear', 0);
        mExtRecvInterp(iPrf, :, iCh) = interp1(currentAlt, mExtRecv(iPrf, 1:length(currentAlt), iCh), altArr, 'linear', 0);
    end

    for iWL = 1:length(lidarCfg.generalCfg.laserWavelengths)
        mExtEmitInterp(iPrf, :, iWL) = interp1(currentAlt, mExtEmit(iPrf, 1:length(currentAlt), iWL), altArr, 'linear', 0);
    end
end
% for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
% 
%     chLabel = lidarCfg.generalCfg.channelLabels{iCh};
% 
%     if NPRF > 1
%         mBscInterp(:, :, iCh) = interp2(meteorDataClean.time, meteorDataClean.altitude, sum(filterFunc(mBscWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM) .* mBsc, 3), mTime, altArr, 'linear');
%         mExtRecvInterp(:, :, iCh) = interp2(meteorDataClean.time, meteorDataClean.altitude, mExtRecv(:, :, iCh), mTime, altArr, 'linear');
%     else
%         mBscInterp(:, :, iCh) = repmat(interp1(meteorDataClean.altitude, reshape(sum(filterFunc(squeeze(mBscWL), lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM) .* squeeze(mBsc(1, :, :)), 2), 1, []), altArr, 'linear'), length(mTime), 1);
%         mExtRecvInterp(:, :, iCh) = repmat(interp1(meteorDataClean.altitude, squeeze(mExtRecv(1, :, iCh)), altArr, 'linear'), length(mTime), 1);
%     end
% 
% end
% 
% for iWL = 1:length(lidarCfg.generalCfg.laserWavelengths)
%     if NPRF > 1
%         mExtEmitInterp(:, :, iWL) = interp2(meteorDataClean.time, meteorDataClean.altitude, mExtEmit(:, :, iWL), mTime, altArr, 'linear');
%     else
%         mExtEmitInterp(:, :, iWL) = repmat(interp1(meteorDataClean.altitude, squeeze(mExtEmit(1, :, iWL)), altArr), length(mTime), 1);
%     end
% end

mBsc = [];
mdr = 0.004;
for iCh = 1:length(lidarCfg.generalCfg.channelLabels)

    chLabel = lidarCfg.generalCfg.channelLabels{iCh};

    switch lidarCfg.generalCfg.channelCfg.(chLabel).channelType
    case 'p'
        mBsc = cat(3, mBsc, mBscInterp(:, :, iCh) * 1 / (1 + mdr));
    case 's'
        mBsc = cat(3, mBsc, mBscInterp(:, :, iCh) * mdr / (1 + mdr));
    case 'e'
        mBsc = cat(3, mBsc, mBscInterp(:, :, iCh));
    case 'r'
        mBsc = cat(3, mBsc, mBscInterp(:, :, iCh));
    otherwise
        error('Unknown channel type: %s', lidarCfg.generalCfg.channelCfg.(chLabel).channelType);
    end

end

%% Calculate Solar Background
solarOpt = getSolarOpt([obsCfg.geoCfg.longitude, obsCfg.geoCfg.latitude], mTime, lidarCfg.generalCfg.zenithAngle, 'spectrumLimits', [300, 1500], 'model', obsCfg.skyCfg.model);

irradiance = NaN(length(mTime), length(lidarCfg.generalCfg.channelLabels));
[solarOptWL, solarOptTIME] = meshgrid(solarOpt.wavelength, solarOpt.time);
for iCh = 1:length(lidarCfg.generalCfg.channelLabels)

    chLabel = lidarCfg.generalCfg.channelLabels{iCh};

    % interpolate solar spectrum
    solarWLFine = solarOpt.wavelength(1):((lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM) / 20):solarOpt.wavelength(end);
    if length(solarOpt.time) > 1
        [TIME, WLFINE] = meshgrid(solarWLFine, mTime);
        solarIrradFine = interp2(solarOptWL, solarOptTIME, solarOpt.irradiance, TIME, WLFINE, 'linear');
    else
        solarIrradFine = reshape(interp1(solarOpt.wavelength, solarOpt.irradiance, solarWLFine, 'linear'), length(mTime), length(solarWLFine));
    end

    % calculate solar background for lidar receiving channel
    switch lidarCfg.generalCfg.channelCfg.(chLabel).channelType
    case 'p'
        irradiance(:, iCh) = trapz(solarWLFine, repmat(filterFunc(solarWLFine, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM), length(mTime), 1) .* solarIrradFine, 2) / 2;
    case 's'
        irradiance(:, iCh) = trapz(solarWLFine, repmat(filterFunc(solarWLFine, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM), length(mTime), 1) .* solarIrradFine, 2) / 2;
    case 'e'
        irradiance(:, iCh) = trapz(solarWLFine, repmat(filterFunc(solarWLFine, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM), length(mTime), 1) .* solarIrradFine, 2);
    case 'r'
        irradiance(:, iCh) = trapz(solarWLFine, repmat(filterFunc(solarWLFine, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM), length(mTime), 1) .* solarIrradFine, 2);
    otherwise
        error('Unknown channel type: %s', lidarCfg.generalCfg.channelCfg.(chLabel).channelType);
    end
end

%% Read Aerosol Scattering Database
emitWLs = [];
for iWL = 1:length(lidarCfg.generalCfg.laserWavelengths)
    emitWLs = cat(2, emitWLs, lidarCfg.generalCfg.laserWavelengths{iWL});
end
recvWLs = [];
for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
    chLabel = lidarCfg.generalCfg.channelLabels{iCh};
    recvWLs = cat(2, recvWLs, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL);
end
aerEmitOpt = getAerOpt(mTime, heightArr, emitWLs, 'aerOptModel', obsCfg.aerCfg.scheme, 'aerOptFilepath', obsCfg.aerCfg.aerOptFilepath);
aerRecvOpt = getAerOpt(mTime, heightArr, recvWLs, 'aerOptModel', obsCfg.aerCfg.scheme, 'aerOptFilepath', obsCfg.aerCfg.aerOptFilepath);

aBsc = [];
for iCh = 1:length(lidarCfg.generalCfg.channelLabels)

    chLabel = lidarCfg.generalCfg.channelLabels{iCh};
    exciteWLInd = lidarCfg.generalCfg.channelCfg.(chLabel).exciteWLInd;

    switch lidarCfg.generalCfg.channelCfg.(chLabel).channelType
    case 'p'
        % Parallel channel
        aBsc = cat(3, aBsc, filterFunc(lidarCfg.generalCfg.laserWavelengths{exciteWLInd}, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM) .* aerRecvOpt.backscatter(:, :, iCh) .* 1 ./ (1 + aerRecvOpt.depolarization(:, :, iCh)));
    case 's'
        % Cross channel
        aBsc = cat(3, aBsc, filterFunc(lidarCfg.generalCfg.laserWavelengths{exciteWLInd}, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM) .* aerRecvOpt.backscatter(:, :, iCh) .* aerRecvOpt.depolarization(:, :, iCh) ./ (1 + aerRecvOpt.depolarization(:, :, iCh)));
    case 'e'
        % Elastic channel
        aBsc = cat(3, aBsc, filterFunc(lidarCfg.generalCfg.laserWavelengths{exciteWLInd}, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM) .* aerRecvOpt.backscatter(:, :, iCh));
    case 'r'
        % Raman channel
        aBsc = cat(3, aBsc, filterFunc(lidarCfg.generalCfg.laserWavelengths{exciteWLInd}, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM) .* aerRecvOpt.backscatter(:, :, iCh));
    otherwise
        error('Unknown channel type: %s', lidarCfg.generalCfg.channelCfg.(chLabel).channelType);
    end

end

%% Calculate Overlap Function
olFactor = NaN(length(lidarCfg.generalCfg.laserWavelengths), length(distArr));
olHeight = NaN(1, length(lidarCfg.generalCfg.laserWavelengths));
for iWL = 1:length(lidarCfg.generalCfg.laserWavelengths)
    teleDiv = (lidarCfg.generalCfg.apertureDiameter / lidarCfg.generalCfg.teleFL);
    olFactor(iWL, :) = overlapCalc(distArr, lidarCfg.generalCfg.laserDivergence{iWL} * 1e-3, teleDiv, lidarCfg.generalCfg.teleDiameter / 2 * 1e-3, lidarCfg.generalCfg.teleLaserDist{iWL} * 1e-3, 0, 0, 0);

    idx = find(olFactor(iWL, :) > 0.96, 1);
    if ~ isempty(idx)
        olHeight(iWL) = heightArr(idx);
    else
        olHeight(iWL) = 0;
    end
end

%% Lidar Signal Simulation
lidarSig = struct();
lidarSig.level0 = struct();   % level 0 signal (with system effects)
lidarSig.level0.time = mTime;
lidarSig.level0.sig = [];
lidarSig.level1 = struct();   % level 1 signal (with no system effects)
lidarSig.level1.altitude = [];
lidarSig.level1.time = mTime;
lidarSig.level1.height = [];
lidarSig.level1.range = [];
lidarSig.level1.sig = [];
lidarSig.level1.bg = [];
trueSig = [];
NtotPoiss = NaN(length(mTime), length(heightArr));
NrawPoiss = NaN(length(mTime), lidarCfg.generalCfg.binNum);
NrawWithPretrigger = NrawPoiss;

for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
    lidarSig.level1.height = heightArr;
    lidarSig.level1.altitude = heightArr + obsCfg.geoCfg.altitude;
    lidarSig.level1.range = distArr;

    chLabel = lidarCfg.generalCfg.channelLabels{iCh};
    exciteWLInd = lidarCfg.generalCfg.channelCfg.(chLabel).exciteWLInd;

    %% Measurement Type
    switch obsCfg.measCfg.type
    case 'normal'

        thisMBsc = mBsc(:, :, iCh);
        thisABsc = aBsc(:, :, iCh);

    case 'depol-cali'

        switch lidarCfg.generalCfg.channelCfg.(chLabel).channelType
        case 'p'
            thisMBsc = mBscInterp(:, :, iCh) / 2;
            thisABsc = filterFunc(lidarCfg.generalCfg.laserWavelengths{exciteWLInd}, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM) .* aerRecvOpt.backscatter(:, :, iCh) / 2;
        case 's'
            thisMBsc = mBscInterp(:, :, iCh) / 2;
            thisABsc = filterFunc(lidarCfg.generalCfg.laserWavelengths{exciteWLInd}, lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL, lidarCfg.generalCfg.channelCfg.(chLabel).filterFWHM) .* aerRecvOpt.backscatter(:, :, iCh) / 2;
        case 'e'
            thisMBsc = mBsc(:, :, iCh);
            thisABsc = aBsc(:, :, iCh);
        case 'r'
            thisMBsc = mBsc(:, :, iCh);
            thisABsc = aBsc(:, :, iCh);
        otherwise
            error('Unknown channel type: %s', lidarCfg.generalCfg.channelCfg.(chLabel).channelType);
        end

    otherwise
        error('Unknown measurement type: %s', obsCfg.measCfg);
    end

    Ns = repmat(olFactor(exciteWLInd, :), length(mTime), 1) .* ...
        lidarCfg.generalCfg.channelCfg.(chLabel).filterPeakTrans .* ...
        lidarCfg.generalCfg.channelCfg.(chLabel).optEffi .* ...
        lidarCfg.generalCfg.channelCfg.(chLabel).quantumEffi * ...
        lidarCfg.generalCfg.laserEnergy{exciteWLInd} * 1e-3 * const.c * ...
        pi * (lidarCfg.generalCfg.teleDiameter / 2 * 1e-3)^2 .* ...
        (thisMBsc + thisABsc) .* ...
        exp(- cumsum((aerEmitOpt.extinction(:, :, exciteWLInd) + mExtEmitInterp(:, :, exciteWLInd) + aerRecvOpt.extinction(:, :, iCh) + mExtRecvInterp(:, :, iCh)) .* repmat([distArr(1), diff(distArr)], length(mTime), 1), 2, 'omitnan')) ./ ...
        repmat(distArr, length(mTime), 1).^2 / 2 * lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL * 1e-9 / (const.c * const.h);
    Nb = irradiance(iCh) * ...
        lidarCfg.generalCfg.channelCfg.(chLabel).filterPeakTrans .* ...
        lidarCfg.generalCfg.channelCfg.(chLabel).optEffi .* ...
        lidarCfg.generalCfg.channelCfg.(chLabel).quantumEffi * ...
        pi * (lidarCfg.generalCfg.apertureDiameter / lidarCfg.generalCfg.teleFL / 2) ^ 2 * ...
        pi * (lidarCfg.generalCfg.teleDiameter * 1e-3 / 2) ^ 2 * ...
        lidarCfg.generalCfg.channelCfg.(chLabel).filterCWL * 1e-9 / (const.c * const.h);
    Nd = lidarCfg.generalCfg.channelCfg.(chLabel).darkCount;

    Ntot = (Ns + Nb + Nd) * lidarCfg.generalCfg.accShots * lidarCfg.generalCfg.binWidth * 1e-9;
    bg = (Nb + Nd) * lidarCfg.generalCfg.accShots * lidarCfg.generalCfg.binWidth * 1e-9;

    for iPrf = 1:length(mTime)
        NtotPoiss(iPrf, :) = sigGenWithNoise(Ntot(iPrf, :), sqrt(Ntot(iPrf, :)), 1, 'poisson');
    end

    lidarSig.level1.sig = cat(3, lidarSig.level1.sig, NtotPoiss);
    lidarSig.level1.bg = cat(3, lidarSig.level1.bg, bg);

    trueSig = cat(3, trueSig, Ntot - bg);

    Nraw = (Ns + Nb + Nd) ./ (1 + lidarCfg.generalCfg.channelCfg.(chLabel).deadtime * 1e-9 * (Ns + Nb + Nd));
    NrawWithPretrigger(:, 1:lidarCfg.generalCfg.nPretrigger) = (Nb + Nd) ./ (1 + lidarCfg.generalCfg.channelCfg.(chLabel).deadtime * 1e-9 * (Nb + Nd));
    NrawWithPretrigger(:, (lidarCfg.generalCfg.nPretrigger + 1):end) = Nraw;

    for iPrf = 1:length(mTime)
        NrawPoiss(iPrf, :) = sigGenWithNoise(NrawWithPretrigger(iPrf, :), sqrt(NrawWithPretrigger(iPrf, :)), 1, 'poisson');
    end

    lidarSig.level0.sig = cat(3, lidarSig.level0.sig, NrawPoiss);

end

%% Output
if ~ isempty(p.Results.savePath)

    if ~ exist(p.Results.savePath, 'dir')
        fprintf('Create save path forcefully.!\n');
        mkdir(p.Results.savePath);
    end

    % raw signal; signal without noise; aerosol and molecular optical properties; overlap factor
    [~, lidarCfgFilename, ext] = fileparts(lidarCfgFile);
    dstFile = fullfile(p.Results.savePath, [lidarCfgFilename, ext]);
    if (exist(dstFile, 'file') ~= 2) && (exist(lidarCfgFile, 'file') == 2) 
        copyfile(lidarCfgFile, dstFile);
    end

    [~, obsCfgFilename, ext] = fileparts(obsCfgFile);
    dstFile = fullfile(p.Results.savePath, [obsCfgFilename, ext]);
    if (exist(dstFile, 'file') ~= 2) && (exist(obsCfgFile, 'file') == 2)
        copyfile(obsCfgFile, dstFile);
    end

    % simulated signal
    rawSigPath = fullfile(p.Results.savePath, 'raw-signal-ASCII');
    if ~ exist(rawSigPath, 'dir')
        mkdir(rawSigPath);
    end

    rawSigBINPath = fullfile(p.Results.savePath, 'raw-signal_BIN');
    if ~ exist(rawSigBINPath, 'dir')
        mkdir(rawSigBINPath);
    end

    for iPrf = 1:length(mTime)

        % write BINARY files
        if length(lidarCfg.generalCfg.laserWavelengths) >= 1
            wavelength1 = lidarCfg.generalCfg.laserWavelengths{1};
            wavelength2 = 0;
            wavelength3 = 0;
        elseif length(lidarCfg.generalCfg.laserWavelengths) >= 2
            wavelength1 = lidarCfg.generalCfg.laserWavelengths{1};
            wavelength2 = lidarCfg.generalCfg.laserWavelengths{2};
            wavelength3 = 0;
        elseif length(lidarCfg.generalCfg.laserWavelengths) >= 3
            wavelength1 = lidarCfg.generalCfg.laserWavelengths{1};
            wavelength2 = lidarCfg.generalCfg.laserWavelengths{2};
            wavelength3 = lidarCfg.generalCfg.laserWavelengths{3};
        else
            wavelength1 = 0;
            wavelength2 = 0;
            wavelength3 = 0;
        end

        recWavelength = [];
        recType = [];
        olHeightAll = [];
        for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
            channelLabel = lidarCfg.generalCfg.channelLabels{iCh};

            recWavelength = cat(2, recWavelength, lidarCfg.generalCfg.channelCfg.(channelLabel).filterCWL);

            switch lidarCfg.generalCfg.channelCfg.(channelLabel).channelType
            case 'p'
                recType = cat(2, recType, 1);
            case 's'
                recType = cat(2, recType, 2);
            case 'r'
                recType = cat(2, recType, 3);
            case 'e'
                recType = cat(2, recType, 0);
            otherwise
                error('Unknown channel type: %s', lidarCfg.generalCfg.channelCfg.(channelLabel).channelType);
            end

            olHeightAll = cat(2, olHeightAll, olHeight(lidarCfg.generalCfg.channelCfg.(channelLabel).exciteWLInd));
        end

        write_CMA_L0(rawSigBINPath, '50000', 'ZPY01', 1000, ...
            obsCfg.geoCfg.longitude, ...
            obsCfg.geoCfg.latitude, ...
            obsCfg.geoCfg.altitude, ...
            mTime(iPrf), ...
            mTime(iPrf) + datenum(0, 1, 0, 0, 0, 1 / lidarCfg.generalCfg.laserRepRate * lidarCfg.generalCfg.accShots), ...
            90 - lidarCfg.generalCfg.zenithAngle, ...
            wavelength1, wavelength2, wavelength3, ...
            length(lidarCfg.generalCfg.channelLabels), ...
            ones(1, length(lidarCfg.generalCfg.channelLabels)), ...
            round(recWavelength), recType, ...
            ones(1, length(lidarCfg.generalCfg.channelLabels)) * lidarCfg.generalCfg.binWidth / 200 * 30, ...
            olHeightAll, ...
            lidarCfg.generalCfg.binNum * ones(1, length(lidarCfg.generalCfg.channelLabels)), ...
            transpose(squeeze(lidarSig.level0.sig(iPrf, :, :))));

        % write ASCII files
        fid = fopen(fullfile(rawSigPath, sprintf('raw-signal-%s.txt', datestr(mTime(iPrf), 'yyyymmdd-HHMMSS'))), 'w');

        metaData = struct();
        metaData.line1 = 'channels:';
        metaData.line2 = sprintf('nPretrigger: %d', lidarCfg.generalCfg.nPretrigger);
        metaData.line3 = sprintf('BinWidth (ns): %d', lidarCfg.generalCfg.binWidth);
        for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
            chLabel = lidarCfg.generalCfg.channelLabels{iCh};
            metaData.line1 = [metaData.line1, ' ', chLabel];
        end
        fprintf(fid, '=======================\n');
        fprintf(fid, 'Measurement Time: %s\n', datestr(mTime(iPrf), 'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid, [metaData.line1, '\n']);
        fprintf(fid, [metaData.line2, '\n']);
        fprintf(fid, [metaData.line3, '\n']);
        fprintf(fid, '=======================\n\n');

        for iBin = 1:size(lidarSig.level0.sig, 2)
            fprintf(fid, [repmat(' %d', 1, size(squeeze(lidarSig.level0.sig(iPrf, :, :)), 2)), '\n'], lidarSig.level0.sig(iPrf, iBin, :));
        end

        fclose(fid);
    end

    % solves
    thisSolves = struct();
    thisSolves.altitude = lidarSig.level1.altitude;
    thisSolves.range = lidarSig.level1.range;
    thisSolves.time = mTime;
    thisSolves.height = lidarSig.level1.height;
    thisSolves.trueSig = trueSig;
    thisSolves.aBsc = aBsc;
    thisSolves.aExt = aerRecvOpt.extinction;
    thisSolves.pdr = aerRecvOpt.depolarization;
    thisSolves.mBsc = mBsc;
    thisSolves.mExt = mExtRecvInterp;
    thisSolves.mdr = mdr;

    [HEIGHT, TIME] = meshgrid(thisSolves.height, thisSolves.time);
    if size(meteorDataClean.time, 1) > 1
        thisSolves.temperature = interp2(meteorDataClean.altitude, meteorDataClean.time, meteorDataClean.temperature, HEIGHT, TIME, 'linear');
        thisSolves.pressure = interp2(meteorDataClean.altitude, meteorDataClean.time, meteorDataClean.pressure, HEIGHT, TIME, 'linear');
        thisSolves.wvmr = interp2(meteorDataClean.altitude, meteorDataClean.time, meteorDataClean.water_vapor_mixing_ratio, HEIGHT, TIME, 'linear');
    else
        thisSolves.temperature = interp1(meteorDataClean.altitude, meteorDataClean.temperature, thisSolves.altitude, 'linear');
        thisSolves.pressure = interp1(meteorDataClean.altitude, meteorDataClean.pressure, thisSolves.altitude, 'linear');
        thisSolves.wvmr = interp1(meteorDataClean.altitude, meteorDataClean.water_vapor_mixing_ratio, thisSolves.altitude, 'linear');
    end

    save(fullfile(p.Results.savePath, 'lidar-solves.mat'), 'thisSolves', 'lidarSig');

end

%% Display

if strcmp(p.Results.visible, 'on')
    %% raw signal for all
    figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', p.Results.visible);

    hold on;

    lineInstances = [];
    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        thisLine = plot(1:lidarCfg.generalCfg.binNum, sum(lidarSig.level0.sig(:, :, iCh), 1), 'LineWidth', 1, 'DisplayName', lidarCfg.generalCfg.channelLabels{iCh});
        lineInstances = cat(2, lineInstances, thisLine);
    end

    hold off;

    xlim([0, lidarCfg.generalCfg.binNum]);
    ylim([1e-1, 1e9]);

    xlabel('Range Bin');
    ylabel('Signal (photon count)');
    title('Lidar Raw Signal');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'FontSize', 11, 'YTick', 10.^(-1:2:9), 'Box', 'on', 'YScale', 'log');

    l = legend(lineInstances, 'Location', 'NorthEast');
    l.FontSize = 10;
    l.Interpreter = 'none';

    if ~ isempty(p.Results.savePath)
        export_fig(gcf, fullfile(p.Results.savePath, 'lidar-raw-signal.png'), '-r300');
    end

    %% effective backscatter for all
    figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', p.Results.visible);

    hold on;

    lineInstances = [];
    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        thisLine = plot(lidarSig.level1.range, sum(lidarSig.level1.sig(:, :, iCh), 1) - sum(lidarSig.level1.bg(:, iCh)), 'LineWidth', 1, 'DisplayName', lidarCfg.generalCfg.channelLabels{iCh});
        lineInstances = cat(2, lineInstances, thisLine);
    end

    hold off;

    xlim([0, max(lidarSig.level1.range)]);
    ylim([1e-1, 1e9]);

    xlabel('Distance (m)');
    ylabel('Signal (photon count)');
    title('Lidar Backscatter Signal');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'YTick', 10.^(-1:2:9), 'FontSize', 11, 'Box', 'on', 'YScale', 'log');

    l = legend(lineInstances, 'Location', 'NorthEast');
    l.FontSize = 10;
    l.Interpreter = 'none';

    if ~ isempty(p.Results.savePath)
        export_fig(gcf, fullfile(p.Results.savePath, 'lidar-backscatter-signal.png'), '-r300');
    end

    %% effective attenuated backscatter for all
    figure('Position', [0, 10, 600, 300], 'Units', 'Pixels', 'Color', 'w', 'visible', p.Results.visible);

    hold on;

    lineInstances = [];
    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        thisLine = plot(lidarSig.level1.range, (sum(lidarSig.level1.sig(:, :, iCh), 1) - sum(lidarSig.level1.bg(:, iCh))) .* lidarSig.level1.range.^2, 'LineWidth', 1, 'DisplayName', lidarCfg.generalCfg.channelLabels{iCh});
        lineInstances = cat(2, lineInstances, thisLine);
    end

    hold off;

    xlim([0, max(lidarSig.level1.range)]);
    ylim([1e8, 1e16]);

    xlabel('Distance (m)');
    ylabel('RCS (a.u.)');
    title('Lidar Range Corrected Signal');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'YTick', 10.^(6:2:13), 'FontSize', 11, 'Box', 'on', 'YScale', 'log');

    l = legend(lineInstances, 'Location', 'NorthEast');
    l.FontSize = 10;
    l.Interpreter = 'none';

    if ~ isempty(p.Results.savePath)
        export_fig(gcf, fullfile(p.Results.savePath, 'lidar-range-corrected-signal.png'), '-r300');
    end

    %% aerosol optical properties
    figure('Position', [0, 10, 600, 600], 'Units', 'Pixels', 'Color', 'w', 'visible', p.Results.visible);

    subpos = subfigPos([0.1, 0.08, 0.87, 0.88], 3, 1, 0, 0.03);

    h = subplot('Position', subpos(1, :), 'Units', 'Normalized');

    hold on;

    lineInstances = [];
    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        thisLine = plot(lidarSig.level1.range, aBsc(1, :, iCh), 'LineWidth', 1, 'DisplayName', lidarCfg.generalCfg.channelLabels{iCh});
        lineInstances = cat(2, lineInstances, thisLine);
    end

    hold off;

    xlim([0, max(lidarSig.level1.range)]);
    ylim([1e-8, 1e-3]);

    xlabel('');
    ylabel('backscatter (m-1sr-1)');
    title('Aerosol Optical Properties');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'YTick', 10.^(-8:2:-3), 'FontSize', 11, 'Box', 'on', 'YScale', 'log', 'XTickLabel', '');

    legendflex(h, lidarCfg.generalCfg.channelLabels, 'ncol', 4, 'interpreter', 'none');

    subplot('Position', subpos(2, :), 'Units', 'Normalized');

    hold on;

    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        plot(lidarSig.level1.range, aerRecvOpt.extinction(1, :, iCh), 'LineWidth', 1, 'DisplayName', lidarCfg.generalCfg.channelLabels{iCh});
    end

    hold off;

    xlim([0, max(lidarSig.level1.range)]);
    ylim([1e-7, 1e-2]);

    xlabel('');
    ylabel('extinction (m-1)');
    title('');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'YTick', 10.^(-7:2:-2), 'FontSize', 11, 'Box', 'on', 'YScale', 'log', 'XTickLabel', '');

    subplot('Position', subpos(3, :), 'Units', 'Normalized');

    hold on;

    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        plot(lidarSig.level1.range, aerRecvOpt.depolarization(1, :, iCh), 'LineWidth', 1, 'DisplayName', lidarCfg.generalCfg.channelLabels{iCh});
    end

    hold off;

    xlim([0, max(lidarSig.level1.range)]);
    ylim([0, 0.5]);

    xlabel('Range (m)');
    ylabel('depolarization');
    title('');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'YTick', 0:0.1:0.5, 'FontSize', 11, 'Box', 'on');

    if ~ isempty(p.Results.savePath)
        export_fig(gcf, fullfile(p.Results.savePath, 'aerosol-optical-properties.png'), '-r300');
    end

    %% molecular optical properties
    figure('Position', [0, 10, 600, 600], 'Units', 'Pixels', 'Color', 'w', 'visible', p.Results.visible);

    subpos = subfigPos([0.1, 0.08, 0.87, 0.88], 3, 1, 0, 0.03);

    subplot('Position', subpos(1, :), 'Units', 'Normalized');

    hold on;

    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        plot(lidarSig.level1.range, mBsc(1, :, iCh), 'LineWidth', 1, 'DisplayName', lidarCfg.generalCfg.channelLabels{iCh});
    end

    hold off;

    xlim([0, max(lidarSig.level1.range)]);
    ylim([1e-14, 1e-4]);

    xlabel('');
    ylabel('backscatter (m-1sr-1)');
    title('Molecular Optical Properties');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'YTick', 10.^(-14:3:-4), 'FontSize', 11, 'Box', 'on', 'YScale', 'log', 'XTickLabel', '');

    subplot('Position', subpos(2, :), 'Units', 'Normalized');

    hold on;

    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        plot(lidarSig.level1.range, mExtRecvInterp(1, :, iCh), 'LineWidth', 1, 'DisplayName', lidarCfg.generalCfg.channelLabels{iCh});
    end

    hold off;

    xlim([0, max(lidarSig.level1.range)]);
    ylim([1e-8, 1e-2]);

    xlabel('');
    ylabel('extinction (m-1)');
    title('');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'YTick', 10.^(-8:2:-2), 'FontSize', 11, 'Box', 'on', 'YScale', 'log', 'XTickLabel', '');

    h = subplot('Position', subpos(3, :), 'Units', 'Normalized');

    hold on;

    lineInstances = [];
    for iCh = 1:length(lidarCfg.generalCfg.channelLabels)
        thisLine = plot(lidarSig.level1.range, mdr * ones(size(lidarSig.level1.range)), 'LineWidth', 1, 'DisplayName', lidarCfg.generalCfg.channelLabels{iCh});
        
        lineInstances = cat(2, lineInstances, thisLine);
    end

    hold off;

    xlim([0, max(lidarSig.level1.range)]);
    ylim([0, 0.1]);

    xlabel('Range (m)');
    ylabel('depolarization');
    title('');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'YTick', 0:0.02:0.1, 'FontSize', 11, 'Box', 'on');

    legendflex(h, lidarCfg.generalCfg.channelLabels, 'ncol', 4, 'interpreter', 'none');

    if ~ isempty(p.Results.savePath)
        export_fig(gcf, fullfile(p.Results.savePath, 'molecule-optical-properties.png'), '-r300');
    end

    %% temperature + pressure + water vapor
    figure('Position', [0, 10, 600, 600], 'Units', 'Pixels', 'Color', 'w', 'visible', p.Results.visible);

    subpos = subfigPos([0.1, 0.08, 0.87, 0.88], 3, 1, 0, 0.03);

    subplot('Position', subpos(1, :), 'Units', 'Normalized');

    hold on;

    plot(meteorDataClean.altitude(1, :), meteorDataClean.temperature(1, :),'Color', 'r');

    hold off;

    xlim([0, max(lidarSig.level1.range)]);
    ylim([-80, 35]);

    xlabel('');
    ylabel('Temperature (\circC)');
    title('Meteorlogical Profiles');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'FontSize', 11, 'Box', 'on', 'XTickLabel', '');

    subplot('Position', subpos(2, :), 'Units', 'Normalized');

    hold on;

    plot(meteorDataClean.altitude(1, :), meteorDataClean.pressure(1, :),'Color', 'g')

    hold off;

    xlim([0, max(lidarSig.level1.range)]);

    xlabel('');
    ylabel('Pressure (hPa)');
    title('');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'FontSize', 11, 'Box', 'on', 'YScale', 'log', 'XTickLabel', '');

    subplot('Position', subpos(3, :), 'Units', 'Normalized');

    hold on;

    plot(meteorDataClean.altitude(1, :), meteorDataClean.water_vapor_mixing_ratio(1, :),'Color', 'b');

    hold off;

    xlim([0, max(lidarSig.level1.range)]);

    xlabel('Range (m)');
    ylabel('Water vapor mixing ratio (g/kg)');
    title('');

    set(gca, 'XMinorTick', 'off', 'YMinorTick', 'on', 'FontSize', 11, 'Box', 'on');

    if ~ isempty(p.Results.savePath)
        export_fig(gcf, fullfile(p.Results.savePath, 'meteorological-profiles.png'), '-r300');
    end
end

end