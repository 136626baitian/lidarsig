global LISAR_ENVS;

LISAR_ENVS.RootDir = fileparts(fullfile(mfilename('fullpath')));
LISAR_ENVS.Authors = 'Zhenping Yin';
LISAR_ENVS.UpdateTime = '2024-01-26';
LISAR_ENVS.Version = '0.0.1';

addpath(genpath(fullfile(LISAR_ENVS.RootDir, 'lib')));
addpath(genpath(fullfile(LISAR_ENVS.RootDir, 'include')));
addpath(genpath(fullfile(LISAR_ENVS.RootDir, 'data')));