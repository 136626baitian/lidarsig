%% 雨滴光学性质模拟 多时间点多高度层计算
clc; clear;
global LISAR_ENVS;

%% 参数设置
savePath = fullfile(LISAR_ENVS.RootDir, 'opac', 'Meteorology_result');
ncfile = 'C:\Users\dong\Desktop\天空背景辐射\cams_shuju\武汉市water2024.11.17\data_stream-oper_stepType-instant.nc';

%% 读取维度信息
mixratio = ncread(ncfile, 'crwc', [1 1 1 1], [1 1 Inf Inf]);  % [lon=1, lat=1, plev, time]

%% 严格保留原始参数定义
water_density = 1000;        
dry_air_mass = 1.29;         
m = 1.33+1e-7i;              
v_avg = 5;                   
wavelengths = [355, 532, 1064];
r = 10.^(linspace(log10(0.05), log10(5), 100)); % 严格保持0.05-5mm定义

%% 预分配结果数组（三维：波长 × 高度层 × 时间）
[nPressure, nTime] = size(mixratio, 3:4); % 获取压力层和时间点数量
ext_results = zeros(length(wavelengths), nPressure, nTime); 

%% 主计算循环（严格保持原始计算步骤）
for iTime = 1:nTime
    for iPressure = 1:nPressure
        % 获取当前混合比（保持原单位转换）
        current_mixratio = mixratio(1, 1, iPressure, iTime);
        lwc = current_mixratio * dry_air_mass;
        I_rain = lwc * v_avg / water_density * 3.6 * 1e6 * 100;
        lambda = 4.1 * I_rain^(-0.21);
        n = 8e3 * exp(-r * lambda);
        
        % 波长循环（保持原始Mie计算）
        for iWL = 1:length(wavelengths)
            ext_rain = zeros(1, length(r));
            k0 = 2 * pi / (wavelengths(iWL)); % 波长单位：纳米
            
            for iR = 1:length(r)
                a = r(iR) * 1000 * 1000; 
                res = Mie(m, k0 * a);     
                ext_rain(iR) = res(1) * pi * (r(iR) * 1e-3)^2;
            end
            
          
            ext_results(iWL, iPressure, iTime) = trapz(r*1e-3, ext_rain .* n);
        end
    end
    fprintf('时间点 %d/%d 完成\n', iTime, nTime);
end

%% 结果保存
save(fullfile(savePath, 'ext_3D_origUnits.mat'), 'ext_results', 'wavelengths', '-v7.3');
disp('计算完成，原始单位版本结果已保存');