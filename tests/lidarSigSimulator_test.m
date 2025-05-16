% test the function of lidar signal simulation
% Zhenping Yin
% 2024-06-14

global LISAR_ENVS;

lidarSigSimulator(datenum(2018, 1, 1, 0, 1, 0), fullfile(LISAR_ENVS.RootDir, 'lib', 'config', 'lidarCfg_default.yml'));
lidarSigSimulator(datenum(2018, 1, 1, 0, 1, 0), fullfile(LISAR_ENVS.RootDir, 'lib', 'config', 'lidarCfg_default.yml'), 'savePath', fullfile(LISAR_ENVS.RootDir, 'tmp'), 'visible', 'on');