%% 云滴光学性质模拟  单时间点测试
clc;
close all;
global LISAR_ENVS;

%% 参数设置
savePath = fullfile(LISAR_ENVS.RootDir, 'opac', 'Meteorology_result');
% 读取NetCDF文件维度信息
ncfile = 'C:\Users\dong\Desktop\天空背景辐射\cams_shuju\武汉市water2024.11.17\data_stream-oper_stepType-instant.nc';

% 读取基础维度数据(可用可不用)
lon = ncread(ncfile, 'longitude');    % 经度
lat = ncread(ncfile, 'latitude');     % 纬度
plev = ncread(ncfile, 'pressure_level');  % 压力层（高度）
unix_time = ncread(ncfile, 'valid_time'); % 读取Unix时间戳（秒）
matlab_time = datetime(unix_time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
mixratio=ncread(ncfile, 'clwc', [1,1,8,1], [1,1,1,1]);%读取混合比(kg/kg)

%赫尔基安-马津公式 云滴粒径谱：𝑛(𝑟)=𝑎𝑟^2*exp⁡(−𝑏𝑟)
water_density=1000; %水的质量密度（1000kg/m^3）
dry_air_mass = 1.29; % 干空气的质量（kg/m^3）
r_q=10;  % 设定r_q平均半径为10微米
b=3/r_q;
m = 1.33+1e-9i; %云滴复折射率
q_w=mixratio*dry_air_mass; 
r = 10.^(linspace(log10(1), log10(20), 100)); % 1-20um
a_param=1.45*1e-6*(q_w/water_density/(r_q^6)); %𝑞_𝑤:单位体积液态水总量 ;𝑞_w=混合比*干空气质量,a命名为a_param，防止重名
wavelengths = [355, 532, 1064];   % incident wavelength (nm)
n=a_param.*r.^2.*exp(-b.*r); %粒径谱
% 计算当前粒径谱对应的q_w
q_w_calc = trapz(r*1e-6, (4/3)*pi*(r*1e-6).^3 * water_density .* n).*1e30;
% 归一化因子
C = q_w / q_w_calc;
% 把n(r)整体归一化
n = C * n;

% 计算消光效率因子
ext_cloud=zeros(1, length(wavelengths));
bsc_cloud=zeros(1, length(wavelengths));
for iWL = 1:length(wavelengths)
    for iR = 1:length(r)
        k0 = 2 * pi / (wavelengths(iWL));
        a = r(iR) * 1000;
        res = Mie(m, k0 * a);
        ext_cloud( iWL, iR) = res(1) * pi * (r(iR) * 1e-6)^2;
        bsc_cloud( iWL, iR) = res(4) * pi * (r(iR) * 1e-6)^2;
    end
end

%%光学性质计算
ext=zeros(1, length(wavelengths));
bsc=zeros(1, length(wavelengths));
for iWL = 1:length(wavelengths)
        ext(iWL) = trapz(r*1e-6 , ext_cloud(iWL, :) .* n.*1e30);
        bsc(iWL) = trapz(r*1e-6 , bsc_cloud(iWL, :) .* n.*1e30)./(4*pi);
end

disp(ext./bsc);
