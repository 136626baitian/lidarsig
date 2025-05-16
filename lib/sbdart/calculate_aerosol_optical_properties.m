function [extinction_coeff, backscatter_coeff] = calculate_aerosol_optical_properties(...
    mixing_ratios, densities, refractive_indices, wavelengths, ...
    particle_sizes, geometric_stds, median_diameters, num_segments)
% 计算气溶胶的消光系数和后向散射系数
% mixing_ratios: 各粒径段的混合比
% densities: 气溶胶密度
% refractive_indices: 各粒径段的折射率
% wavelengths: 入射波长 (nm)
% particle_sizes: 各粒径段的粒径分布 (um)
% geometric_stds: 各粒径段的几何标准差
% median_diameters: 各粒径段的中值直径
% num_segments: 粒径段数量

% 初始化结果
extinction_coeff = zeros(1, length(wavelengths));
backscatter_coeff = zeros(1, length(wavelengths));

% 初始化粒子大小分布
for iSeg = 1:num_segments
    % 定义粒径分布
    seg_ratio = mixing_ratios(iSeg);
    seg_density = densities(iSeg);
    seg_particles = particle_sizes{iSeg};
    seg_std = geometric_stds(iSeg);
    seg_median = median_diameters(iSeg);
    seg_refractive = refractive_indices(iSeg);
    
    % 计算归一化粒子常数
    c = Particle_Constants_Calculation(seg_ratio, seg_density, seg_particles, seg_std, seg_median);
    
    % 计算消光效率因子和后向散射系数
    [ext, bsc] = calculate_segment_optics(c, seg_particles, seg_refractive, wavelengths,seg_std, seg_median);
    
    % 累加到总体结果
    extinction_coeff = extinction_coeff + ext;
    backscatter_coeff = backscatter_coeff + bsc;
end
end

function [ext, bsc] = calculate_segment_optics(c, particles, m, wavelengths,seg_std, seg_median)
% 计算单个粒径段的消光和后向散射
% c: 归一化粒子常数
% particles: 粒径分布
% m: 折射率
% wavelengths: 波长 (nm)

% 初始化结果
ext = zeros(1, length(wavelengths));
bsc = zeros(1, length(wavelengths));

% 计算每个波长的消光和后向散射
for iWL = 1:length(wavelengths)
    k0 = 2 * pi / wavelengths(iWL);
    a = particles * 1000; % 转换为米
    % 对每个粒径进行计算
    ext_i = zeros(size(a));
    bsc_i = zeros(size(a));
    for iPart = 1:length(a)
        x = k0 * a(iPart);
        res = Mie(m, x);
        ext_i(iPart) = res(1) * pi * (a(iPart)* 1e-6)^2; % 消光效率因子
        bsc_i(iPart) = res(4) * pi * (a(iPart)* 1e-6)^2; % 后向散射效率因子
    end
    
    % 计算体积密度谱
    v_sd = (4/3 * pi * (particles).^3) * c / (sqrt(2*pi) * log(seg_std)) .* ...
        exp(-(log(particles) - log(seg_median)).^2 / (2 * (log(seg_std))^2));
    
    % 积分计算总体消光和后向散射
    ext(iWL) = trapz(particles * 1e-6, ext_i .* 1e6 .* v_sd ./ (4/3 * pi * (particles * 1e-6).^3))/1e6;
    bsc(iWL) = trapz(particles * 1e-6, bsc_i .* 1e6 .* v_sd ./ (4/3 * pi * (particles * 1e-6).^3)) /1e6/ (4 * pi);
end
end