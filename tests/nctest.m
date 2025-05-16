% NetCDF 文件路径
ncfile = 'C:\Users\dong\Desktop\lidar\lidar shuju\radiosonde_57494_20180101_0000.nc';

% 读取高度数据
altitude = ncread(ncfile, 'altitude');

% 读取水汽混合比数据
water_vapor_mixing_ratio = ncread(ncfile, 'water_vapor_mixing_ratio');

% 自定义水汽混合比数据
custom_water_vapor_mixing_ratio = [16; 15; 14; 13; 11; 10; 9; 8.8; 9; 8.8; ...
                                   8.5; 8.2; 7.9; 7.6; 7.3; 7.000; 6.737; 6.474; 6.210; 5.947; 5.684; 5.421; 5.158; 4.894; 4.631; 4.368; ...
                                   4.105; 3.842; 3.578; 3.315; 3.052; 2.789; 2.526; 2.262; 2.000; 1.8; 1; 0.8; 0.6; 0.4; 0.2; 0.1; ...
                                   0.05; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; ...
                                   NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; ...
                                   NaN; NaN; NaN; NaN; NaN; NaN; NaN];
% 检查维度是否一致
if ~isequal(size(water_vapor_mixing_ratio), size(custom_water_vapor_mixing_ratio))
    error('自定义数据的维度与 NetCDF 文件中的数据不匹配！');
end

% 修改数据
water_vapor_mixing_ratio = custom_water_vapor_mixing_ratio;

% 将修改后的数据写回 NetCDF 文件
ncwrite(ncfile, 'water_vapor_mixing_ratio', water_vapor_mixing_ratio);

disp('NetCDF 文件中的水汽混合比数据已成功更新！');