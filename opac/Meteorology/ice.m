%% 冰晶光学性质模拟  单时间点测试
clc;
close all;
global LISAR_ENVS;

%% 参数设置
savePath = fullfile(LISAR_ENVS.RootDir, 'opac', 'Meteorology_result');
% 读取NetCDF文件维度信息
ncfile = 'C:\Users\dong\Desktop\天空背景辐射\cams_shuju\武汉市2024.1.21\3ac07dd94fac6957f1c4a36343422604.nc';

% 读取基础维度数据(可用可不用)
lon = ncread(ncfile, 'longitude');    % 经度
lat = ncread(ncfile, 'latitude');     % 纬度
plev = ncread(ncfile, 'pressure_level');  % 压力层（高度）
unix_time = ncread(ncfile, 'valid_time'); % 读取Unix时间戳（秒）
matlab_time = datetime(unix_time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
mixratio=ncread(ncfile, 'ciwc', [1,1,13,11], [1,1,1,1]);%读取混合比(kg/kg)


%% 冰晶粒径谱：𝑛(D)=N_total*C*(D/D_M)^μ*exp(−(4+μ)*D/D_M) C为归一化因子 D为粒子直径
% C经过推导后，公式如下：C=(4+μ)^1+μ/(D_M*Γ(μ+1))
ice_density=917;      %冰的质量密度（917kg/m^3）
dry_air_mass = 1.29; % 干空气的质量（kg/m^3）
mu=2;                 %形状参数,目前假设为2（典型冰晶分布值）
D_M=300;             %假设中位直径为300um
C=(4+mu)^(mu+1)/(D_M*1e-6*gamma(mu + 1));
N_total=mixratio*dry_air_mass*(4+mu)^3*gamma(mu + 1)/((pi/6)*ice_density*(D_M*1e-6)^3*gamma(mu + 4));
m = 1.31+1e-8i; %冰晶复折射率
D = 10.^(linspace(log10(50), log10(1000), 100)); % 50-1000um
wavelengths = [355, 532, 1064];   % incident wavelength (nm)
n=N_total*C*(D/D_M).^mu.*exp(-(4+mu)*D/D_M); %粒径谱

% 计算消光效率因子
ext_cloud=zeros(1, length(wavelengths));
bsc_cloud=zeros(1, length(wavelengths));
for iWL = 1:length(wavelengths)
    for iR = 1:length(D)
        k0 = 2 * pi / (wavelengths(iWL));
        a = D(iR) * 1000;
        res = Mie(m, k0 * a);
        ext_cloud( iWL, iR) = res(1) * pi * (D(iR)/2 * 1e-6)^2;
        bsc_cloud( iWL, iR) = res(4) * pi * (D(iR)/2 * 1e-6)^2;
    end
end


%%光学性质计算
ext=zeros(1, length(wavelengths));
bsc=zeros(1, length(wavelengths));
for iWL = 1:length(wavelengths)
        ext(iWL) = trapz(D*1e-6 , ext_cloud(iWL, :) .* n);
        bsc(iWL) = trapz(D*1e-6 , bsc_cloud(iWL, :) .* n)./(4*pi);
end
disp(ext./bsc)