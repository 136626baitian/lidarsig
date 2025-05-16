%% 雨滴光学性质模拟 单时间点测试 - Gamma分布版
clc;
clear;
close all;
global LISAR_ENVS;

%% 参数设置
savePath = fullfile(LISAR_ENVS.RootDir, 'opac', 'Meteorology_result');

% 读取NetCDF文件
ncfile = 'C:\Users\dong\Desktop\天空背景辐射\cams_shuju\武汉市water2024.11.17\data_stream-oper_stepType-instant.nc';

% 读取基础维度数据
lon = ncread(ncfile, 'longitude');    % 经度
lat = ncread(ncfile, 'latitude');     % 纬度
plev = ncread(ncfile, 'pressure_level');  % 压力层
unix_time = ncread(ncfile, 'valid_time'); % Unix时间戳（秒）
matlab_time = datetime(unix_time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

mixratio = ncread(ncfile, 'crwc', [1,1,8,1], [1,1,1,1]); % 混合比(kg/kg)

%% 基础物理参数
water_density = 1000;    % 水的密度 kg/m^3
dry_air_mass = 1.29;     % 干空气密度 kg/m^3
m = 1.33 + 1e-7i;        % 雨滴复折射率
v_avg = 5;               % 雨滴平均下落速度 m/s
lwc = mixratio * dry_air_mass;   % 液态水含量 kg/m^3
r = 10.^(linspace(log10(0.05), log10(7), 120)); % 粒径范围：0.05mm-7mm
wavelengths = [355, 532, 1064]; % 波长 (nm)

%% 雨强计算
I_rain = lwc * v_avg / water_density * 3.6e6;  % mm/h

%% Gamma分布雨滴谱计算
Lambda = 4.1 * I_rain^(-0.21);       % Lambda (1/mm)
% 选择μ
if I_rain < 1
    mu = 0;
elseif I_rain < 10
    mu = 2;
else
    mu = 3;
end
N_0 = 8e2 * I_rain^(0.232).*r.^(-mu);      % N0 (m^-3 mm^-1)
% 雨滴谱
n = N_0.* (r.^mu) .* exp(-Lambda * r);  % m^-3 mm^-1

% 计算当前粒径谱对应的LWC
lwc_calc = trapz(r*1e-3, (4/3)*pi*(r*1e-3).^3 .* water_density .* n) * 100;
C = lwc / lwc_calc;  % 归一化因子

%% Mie散射计算（严格保留 a = r*1000*1000）
ext_rain = zeros(length(wavelengths), length(r));
bsc_rain = zeros(length(wavelengths), length(r));

for iWL = 1:length(wavelengths)
    for iR = 1:length(r)
        k0 = 2 * pi / wavelengths(iWL); % 波数 (1/nm)
        a = r(iR) * 1000 * 1000;         
        res = Mie(m, k0 * a);             % Mie计算
        ext_rain(iWL, iR) = res(1) * pi * (r(iR) * 1e-3)^2; % 消光截面 (m^2)
        bsc_rain(iWL, iR) = res(4) * pi * (r(iR) * 1e-3)^2; % 后向散射截面 (m^2)
    end
end

%% 光学性质积分
ext = zeros(1, length(wavelengths));
bsc = zeros(1, length(wavelengths));

for iWL = 1:length(wavelengths)
    ext(iWL) = trapz(r*1e-3, ext_rain(iWL, :) .* C .* n.*100);   % 消光系数 (1/m)
    bsc(iWL) = trapz(r*1e-3, bsc_rain(iWL, :) .* C .* n.*100) / (4*pi); % 后向散射系数 (1/m/sr)
end



