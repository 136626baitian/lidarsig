%% 冰晶光学性质模拟 多时间点多高度层计算
clc; clear;
global LISAR_ENVS;

%% 参数设置
savePath = fullfile(LISAR_ENVS.RootDir, 'opac', 'Meteorology_result');
ncfile = 'C:\Users\dong\Desktop\天空背景辐射\cams_shuju\武汉市2024.1.21\3ac07dd94fac6957f1c4a36343422604.nc';

%% 读取全量数据（第一经纬度点，所有压力层和时间）
mixratio = ncread(ncfile, 'ciwc', [1 1 1 1], [1 1 Inf Inf]); % [lon, lat, plev, time]

%% 保持原始参数定义
ice_density = 917;        
dry_air_mass = 1.29;      
mu = 2;                   
D_M = 300;                % 保持300μm中位直径
m = 1.31+1e-8i;           
D = 10.^(linspace(log10(50), log10(1000), 100)); % 严格保持50-1000μm定义
wavelengths = [355, 532, 1064]; 

%% 预分配结果数组
[nPressure, nTime] = size(mixratio, 3:4);
nWL = length(wavelengths);
ext_results = zeros(nWL, nPressure, nTime);  % 消光系数矩阵
bsc_results = zeros(nWL, nPressure, nTime);  % 后向散射系数矩阵

%% 主计算循环（保持原始单位体系）
for iTime = 1:nTime
    for iPressure = 1:nPressure
        % 获取当前时空的混合比
        current_mixratio = mixratio(1, 1, iPressure, iTime);
        
        % --- 严格保持原始计算流程 ---
        C = (4+mu)^(mu+1)/(D_M*1e-6*gamma(mu + 1));
        N_total = current_mixratio * dry_air_mass * (4+mu)^3 * gamma(mu + 1) /...
                 ((pi/6)*ice_density*(D_M*1e-6)^3*gamma(mu + 4));
        n = N_total*C*(D/D_M).^mu.*exp(-(4+mu)*D/D_M);
        
        % 波长循环
        for iWL = 1:nWL
            % 预初始化
            ext_cloud = zeros(1, length(D));
            bsc_cloud = zeros(1, length(D));
            k0 = 2 * pi / (wavelengths(iWL)); % 保持原纳米波长单位
            
            % 粒径循环（保持原始单位转换）
            for iD = 1:length(D)
                a = D(iD) * 1000; % 严格保持原转换：μm → nm
                res = Mie(m, k0 * a); % 保持原Mie参数计算
                radius = D(iD)/2 * 1e-6; % 保持原半径计算（μm → m）
                
                ext_cloud(iD) = res(1) * pi * radius^2;
                bsc_cloud(iD) = res(4) * pi * radius^2;
            end
            
            % 积分计算（保持原方法）
            ext_results(iWL, iPressure, iTime) = trapz(D*1e-6, ext_cloud .* n);
            bsc_results(iWL, iPressure, iTime) = trapz(D*1e-6, bsc_cloud .* n)/(4*pi);
        end
    end
    fprintf('已完成时间点 %d/%d\n', iTime, nTime);
end

%% 结果保存（新增高度和时间维度）
save(fullfile(savePath, 'ice_properties_3D.mat'),...
    'ext_results', 'bsc_results', 'plev', 'matlab_time', 'wavelengths', '-v7.3');
disp('冰晶三维光学性质计算完成');