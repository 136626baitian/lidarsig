function c = Particle_Constants_Calculation(mixing_ratio, aerosol_density, r,sgm,rgm) %输入变量：混合比，气溶胶密度，粒径，几何标准差，峰值半径
    % 常量定义  
    dry_air_mass = 1.29; % 干空气的质量（kg/m^3）  
    
    % 计算单个粒子的质量（kg）
    single_particle_mass = (4/3) * pi * (aerosol_density * (r.^3));

    % 计算单位体积气溶胶的质量（kg/m^3）
    aerosol_mass_per_volume = mixing_ratio * dry_air_mass;

    % 计算归一化粒子常数计算
    c = aerosol_mass_per_volume / aerosol_density * sqrt(2 * pi) / (4 / 3 * pi*trapz(r ,...
        (r.^3)./log(sgm).* exp(- (log(r) - log(rgm)).^2 / (2 * (log(sgm)).^2))));  %c(个/um3)
end