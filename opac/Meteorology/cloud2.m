%% 云滴光学性质模拟  多时间点测试
clc;
clear;
close all;
global LISAR_ENVS;

%% 参数设置
% 文件路径
ncfile = 'C:\Users\dong\Desktop\天空背景辐射\cams_shuju\武汉市water2024.11.17\data_stream-oper_stepType-instant.nc';
savePath = fullfile(pwd, 'optical_properties.nc'); % 结果保存路径

% 物理参数
water_density = 1000;     % 水的密度 [kg/m³]
dry_air_mass = 1.29;      % 干空气质量 [kg/m³]
r_q = 10;                 % 平均半径 [微米]
wavelengths = [355, 532, 1064]; % 波长 [nm]
m = 1.33 + 1e-9i;         % 水滴复折射率

% 粒径分布参数
r = 10.^(linspace(log10(1), log10(20), 100)); % 1-20微米（100个对数间隔点）
b = 3 / r_q;              % 赫尔基安-马津公式参数

%% 读取NetCDF数据
% 基础维度
lon = ncread(ncfile, 'longitude');
lat = ncread(ncfile, 'latitude');
plev = ncread(ncfile, 'pressure_level'); % 压力层 [Pa]
unix_time = ncread(ncfile, 'valid_time');

% 固定经纬度为第一个点
lon_idx = 1;
lat_idx = 1;

% 获取数据维度信息
n_times = length(unix_time);
n_plev = length(plev);

% 读取液态水混合比（维度顺序：lon×lat×height×time）
clwc = ncread(ncfile, 'clwc', [lon_idx, lat_idx, 1, 1], [1, 1, n_plev, n_times]); % [kg/kg]

%% 预计算Mie散射系数
% 初始化效率因子矩阵
Qext = zeros(length(wavelengths), length(r));
Qbsc = zeros(length(wavelengths), length(r));

for iWL = 1:length(wavelengths)
    for iR = 1:length(r)
        k0 = 2 * pi / wavelengths(iWL) ; 
        a = r(iR) * 1000;                        
        res = Mie(m, k0*a);                       % 调用Mie函数
        Qext(iWL, iR) = res(1);                   % 消光效率因子
        Qbsc(iWL, iR) = res(4);                   % 后向散射效率因子
    end
end

%% 主计算循环
ext_all = zeros(n_plev, n_times, length(wavelengths)); % 维度: 高度×时间×波长
bsc_all = zeros(n_plev, n_times, length(wavelengths)); % 维度: 高度×时间×波长

for iTime = 1:n_times
    for iPlev = 1:n_plev
        % 当前混合比
        mixratio = clwc(1, 1, iPlev, iTime);
        
        % 计算液态水含量 [kg/m³]
        q_w = mixratio * dry_air_mass;
        
        % 赫尔基安-马津参数
        a_param = 1.45e-6 * (q_w / water_density) / (r_q^6);
        
        % 生成粒径分布 [个/um⁻³·um⁻¹]
        n_dist = a_param .* r.^2 .* exp(-b*r); 
        
        % 计算光学系数
        for iWL = 1:length(wavelengths)
            % 消光系数计算
            ext_integrand = pi * (r).^2 .* Qext(iWL, :) .* n_dist*1e24;
            ext_all(iPlev, iTime, iWL) = trapz(r, ext_integrand); % [m⁻¹]
            
            % 后向散射系数计算
            bsc_integrand = pi * (r).^2 .* Qbsc(iWL, :) .* n_dist*1e24;
            bsc_all(iPlev, iTime, iWL) = trapz(r, bsc_integrand) / (4*pi); % [m⁻1·sr⁻1]
        end
    end
end


