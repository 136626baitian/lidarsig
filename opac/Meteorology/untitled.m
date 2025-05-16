%% 雨滴光学性质模拟 - Gamma分布版 - 消光随雨强变化
clc;
clear;
close all;
global LISAR_ENVS;

%% 固定参数
% 读取NetCDF这里就略过，假设基础物理量给定
water_density = 1000;    % 水的密度 kg/m^3
dry_air_mass = 1.29;     % 干空气密度 kg/m^3
m = 1.33 + 1e-7i;        % 雨滴复折射率
v_avg = 5;               % 雨滴平均下落速度 m/s
r = 10.^(linspace(log10(0.05), log10(7), 120)); % 粒径范围：0.05mm-7mm
wavelengths = [355, 532, 1064]; % 波长 (nm)

%% 设置雨强范围
rain_rates = logspace(log10(0.1), log10(100), 20); % 0.1到100 mm/h，共20个点
extinction_all = zeros(length(wavelengths), length(rain_rates));

%% 循环遍历每一个雨强
for iRR = 1:length(rain_rates)
    % 通过雨强反推lwc
    I_rain = rain_rates(iRR); % 当前雨强 mm/h
    lwc = I_rain / (v_avg * 3.6e6) * water_density; % 单位 m^3
    
    % Gamma分布参数更新
    Lambda = 4.1 * I_rain^(-0.21); % (1/mm)
    
    % 选择μ
    if I_rain < 1
        mu = 0;
    elseif I_rain < 10
        mu = 2;
    else
        mu = 3;
    end

    % N0 更新（注意这里不乘r^-mu，不改你的思路）
    N_0 = 8e2 * I_rain^(0.232); % m^-3 mm^-1
    
    % 雨滴谱
    n = N_0 * (r.^mu) .* exp(-Lambda * r); % m^-3 mm^-1

    % 归一化（按照液态水含量校准）
    lwc_calc = trapz(r*1e-3, (4/3)*pi*(r*1e-3).^3 .* water_density .* n) * 100; % g/m^3
    C = lwc / lwc_calc;

    % Mie散射计算（保留你的a = r*1000*1000）
    ext_rain = zeros(length(wavelengths), length(r));
    for iWL = 1:length(wavelengths)
        for iR = 1:length(r)
            k0 = 2 * pi / wavelengths(iWL); % 波数 (1/nm)
            a = r(iR) * 1000 * 1000;         
            res = Mie(m, k0 * a);
            ext_rain(iWL, iR) = res(1) * pi * (r(iR) * 1e-3)^2; % m^2
        end
    end

    % 光学积分
    for iWL = 1:length(wavelengths)
        extinction_all(iWL, iRR) = trapz(r*1e-3, ext_rain(iWL, :) .* C .* n.*100); % m^-1
    end
end

%% 绘图
figure;
hold on;
colors = lines(length(wavelengths));
for iWL = 1:length(wavelengths)
    plot(rain_rates, extinction_all(iWL, :), '-o', 'Color', colors(iWL,:), 'DisplayName', sprintf('%d nm', wavelengths(iWL)));
end
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Rain Rate (mm/h)');
ylabel('Extinction Coefficient (m^{-1})');
title('Extinction Coefficient vs Rain Rate (Gamma distribution, constant a scheme)');
legend('show');
grid on;
