% test reading aerosol optical properties from file
% Zhenping Yin
% 2024-08-08

global LISAR_ENVS;

workPath = fileparts(mfilename('fullpath'));

%% Define aerosol optical properties
mTime = datenum(2018, 1, 1, 0, 1, 0);
height = 0:30:40000;
wavelength = [352, 1064];
savePath = 'C:\Users\dong\Desktop\文章材料\压缩包\LiSPP-master\LiSPP-master\analysis\wv-t-retrievals-under-dense-aerosols\data\shiyan2';

% 设置消光系数曲线，
extinction = zeros(1, length(height), length(wavelength));
extinction = reshape(transpose([1e-4 * ones(size(height)); 1e-5 * ones(size(height))]), 1, length(height), length(wavelength));

extinction(1, :, 1) = 0;
extinction(1, :, 2) = 0;
extinction(1, 1:167, 1) = 5e-4;
extinction(1, 1:167, 2) = 8e-5;


backscatter = extinction / 50;
depol = reshape(transpose([0.05 * ones(size(height)); 0.05 * ones(size(height))]), 1, length(height), length(wavelength));
saveAerOpt(fullfile(LISAR_ENVS.RootDir, 'data', 'aerosol_optical_properties.nc'), mTime, height, wavelength, extinction, backscatter, depol);

%% Signal Simulation
lidarSigSimulator(datenum(2018, 1, 1, 0, 1, 0), fullfile(LISAR_ENVS.RootDir, 'lib', 'config', 'TWLidarCfg.yml'), ...
    'obsCfgFile', fullfile(workPath, 'config', 'obsCfg_test.yml'), ...
    'visible', 'on', ...
    'savePath', savePath);
