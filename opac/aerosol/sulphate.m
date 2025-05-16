clc;
close all;
global LISAR_ENVS;
workPath = fileparts(mfilename('fullpath'));
savePath = 'C:\Users\dong\Desktop\lidar\LiSAR-master\result';

%% 硫酸气溶胶光学性质计算
% 读取NetCDF文件维度信息
ncfile = 'C:\Users\dong\Desktop\天空背景辐射\cams_shuju\93ea7e416c9ab17035e3beb650a73df8\data_plev.nc';

% 读取基础维度数据
lon = ncread(ncfile, 'longitude');    % 经度
lat = ncread(ncfile, 'latitude');     % 纬度
plev = ncread(ncfile, 'pressure_level');  % 压力层（高度）
unix_time = ncread(ncfile, 'valid_time'); % 读取Unix时间戳（秒）
matlab_time = datetime(unix_time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

% 获取各维度长度
n_lon = length(lon);
n_lat = length(lat);
n_plev = length(plev);
n_time = length(unix_time);
wavelengths = [355, 532, 1064]; % 波长参数
n_wl = length(wavelengths);

%% 初始化数组（维度顺序：时间×高度×波长）注意时间维度，多个时间需要三维数组
ext_array = zeros(n_plev, n_wl); % 消光系数数组
bsc_array = zeros(n_plev, n_wl); % 后向散射系数数组

%% 亲水性有机物气溶胶固定参数(不同气溶胶需要更改参数)
num_segments = 1;
densities = [1800]; % kg/m³
refractive_indices = [1.48 + 0.001i]; %复折射率
geometric_stds = [1.8];  %几何标准差
median_diameters = [0.3]; % 峰值半径

% 粒径分布定义
r1 = 10.^(linspace(log10(0.03), log10(0.5), 100)); 
particle_sizes = {r1};

%% 主计算循环（仅处理第一个经纬度点）
for i_plev = 1:n_plev
     i_time = 5;    %这里计算12点气溶胶光学性质
        try
            % 固定经纬度索引为1
            i_lon = 1;
            i_lat = 1;

            % 读取混合比（索引固定为第一个经纬度）
            aermr11 = ncread(ncfile, 'aermr11', [i_lon, i_lat, i_plev, i_time], [1,1,1,1]);
            % 计算光学参数
            mixing_ratios = [aermr11];
            [ext, bsc] = calculate_aerosol_optical_properties(...
                mixing_ratios, densities, refractive_indices, wavelengths, ...
                particle_sizes, geometric_stds, median_diameters, num_segments);

            % 存储到二维数组（高度×波长）
            ext_array(i_plev, :) = ext; 
            bsc_array(i_plev, :) = bsc;
        catch ME
            % 错误处理：记录错误并继续执行
            warning('Error at (plev=%d, time=%d): %s', i_plev, i_time, ME.message);
            ext_array(i_time, i_plev, :) = NaN; % 填充NaN
            bsc_array(i_time, i_plev, :) = NaN;
        end
end

disp(ext_array./bsc_array);