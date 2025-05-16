function [sza, cos_sza] = calculate_sza(lat, lon, year, month, day, local_hour, timezone)
% 计算太阳天顶角（SZA）及其余弦值
% 输入参数：
%   lat, lon  - 地理坐标（纬度：北纬为正；经度：东经为正）
%   year      - 本地年份（如2023）
%   month     - 本地月份（1-12）
%   day       - 本地日期（1-31）
%   local_hour- 本地时间的小时（0-24）
%   timezone  - 时区（如北京时间+8）
% 输出参数：
%   sza      - 太阳天顶角（单位：度）
%   cos_sza  - 天顶角余弦值（范围[0,1]）

% 参数校验
validateattributes(lat, {'numeric'}, {'scalar', '>=', -90, '<=', 90});
validateattributes(lon, {'numeric'}, {'scalar', '>=', -180, '<=', 180});
validateattributes(year, {'numeric'}, {'integer', 'scalar', '>=', 1800});
validateattributes(month, {'numeric'}, {'integer', 'scalar', '>=', 1, '<=', 12});
validateattributes(day, {'numeric'}, {'integer', 'scalar', '>=', 1, '<=', 31});
validateattributes(local_hour, {'numeric'}, {'scalar', '>=', 0, '<=', 24});
validateattributes(timezone, {'numeric'}, {'integer', 'scalar', '>=', -12, '<=', 14});

% --- 步骤1: 本地时间转UTC时间 ---
utc_hour = local_hour - timezone;
day_offset = 0;

% 处理小时溢出
if utc_hour < 0
    utc_hour = utc_hour + 24;
    day_offset = -1;
elseif utc_hour >= 24
    utc_hour = utc_hour - 24;
    day_offset = 1;
end

% 初步UTC日期
utc_year = year;
utc_month = month;
utc_day = day + day_offset;

% 处理日期溢出（跨月/年）
if utc_day < 1
    utc_month = utc_month - 1;
    if utc_month < 1
        utc_month = 12;
        utc_year = utc_year - 1;
    end
    utc_day = eomday(utc_year, utc_month); % 获取上月最后一天
elseif utc_day > eomday(utc_year, utc_month)
    utc_day = 1;
    utc_month = utc_month + 1;
    if utc_month > 12
        utc_month = 1;
        utc_year = utc_year + 1;
    end
end

% --- 步骤2: 计算儒略日（修正1月/2月问题）---
% 调整月份为3-14月（对应上一年）
if utc_month <= 2
    calc_year = utc_year - 1;
    calc_month = utc_month + 12;
else
    calc_year = utc_year;
    calc_month = utc_month;
end

a = floor((14 - calc_month)/12);
y = calc_year + 4800 - a;
m = calc_month + 12*a - 3;

JDN = utc_day + floor((153*m + 2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045;
JD = JDN + (utc_hour - 12)/24; % 儒略日含小数

% --- 步骤3: 计算太阳位置 ---
n = JD - 2451545.0;            % 自J2000起的天数
L = mod(280.460 + 0.9856474*n, 360);  % 平黄经
g = mod(357.528 + 0.9856003*n, 360);  % 平近点角
lambda = L + 1.915*sind(g) + 0.020*sind(2*g); % 真黄经

% 赤纬计算
epsilon = 23.439 - 0.0000004*n;        % 黄赤交角
decl = asind(sind(epsilon).*sind(lambda));

% 时角计算
GMST = mod(280.46061837 + 360.98564736629*(JD - 2451545.0) + ...
       0.000387933*( (JD - 2451545.0)/36525 )^2 - ...
       ( (JD - 2451545.0)/36525 )^3 / 38710000, 360);
LST = GMST + lon;              % 本地恒星时
hour_angle = LST - lambda;     % 时角
hour_angle = mod(hour_angle + 180, 360) - 180; % 归一化到[-180,180]

% --- 步骤4: 天顶角计算 ---
lat_rad = deg2rad(lat);
decl_rad = deg2rad(decl);
ha_rad = deg2rad(hour_angle);

cos_sza = sin(lat_rad).*sin(decl_rad) + cos(lat_rad).*cos(decl_rad).*cos(ha_rad);
cos_sza = max(min(cos_sza, 1), 0); % 限制物理范围
sza = acosd(cos_sza);

end