function [solarOpt] = getSolarOpt(geoLoc, measTime, zenithAngle, varargin)
% GETSOLAROPT calculate solar irradiance via SBDART or other models.
%
% USAGE:
%    [output] = getSolarOpt(geoLoc, measTime, zenithAngle)
%
% INPUTS:
%    geoLoc: numeric
%        [longitude (-180, 180), latitude (-90, 90)] in degree
%    measTime: numeric
%        measurement time. (UTC+0)
%    zenithAngle: numeric
%        zenith angle. (degree)
%
% KEYWORDS:
%    spectrumLimits: numeric
%        spectrum limits. (nm)
%    model: char
%        none
%        sbdart
%        daytime
%
% OUTPUTS:
%    solarOpt: struct
%        wavelength (nm)
%        sky spectral irradiance (W m^-2 nm^-1 sr^-1)
%
% EXAMPLE:
%
% HISTORY:
%    2024-03-19: first edition by Zhenping
%    2025-03-8: integrated SBDART model support

global LISAR_ENVS;

% 输入参数解析
p = inputParser;
addRequired(p, 'geoLoc', @isnumeric);          % 经纬度 [经度, 纬度]
addRequired(p, 'measTime', @isnumeric);        % 时间数组
addRequired(p, 'zenithAngle', @isnumeric);     % 天顶角数组（度）
addParameter(p, 'spectrumLimits', [300, 1200], @isnumeric);  % 波长范围（nm）
addParameter(p, 'model', 'none', @ischar);      % 模型选择
parse(p, geoLoc, measTime, zenithAngle, varargin{:});

% 初始化输出结构
solarOpt = struct();
solarOpt.time = measTime(:);  % 确保列向量
solarOpt.wavelength = [];
solarOpt.irradiance = [];

% 模型分支
switch lower(p.Results.model)
    case 'none'
        solarOpt.wavelength = p.Results.spectrumLimits(1):p.Results.spectrumLimits(2);
        solarOpt.irradiance = 0.01 * ones(size(solarOpt.wavelength));

    case 'sbdart'
        %% SBDART模型处理
        solarOpt.wavelength = (p.Results.spectrumLimits(1) : 5 : p.Results.spectrumLimits(2)); % 步长5nm
        solarOpt.irradiance = zeros(length(measTime), length(solarOpt.wavelength)); 
        % 路径配置
        sbdartPath='D:\cygwin\bin\bash.exe';
        sbdartInputDir='C:\Users\dong\Desktop\SBDART-master';
        cygwinSbdartPath = '/cygdrive/c/Users/dong/Desktop/SBDART-master/sbdart.exe';  % 这三个路径需要根据个人电脑路径进行修改
        %% 动态计算天顶角
        lat = geoLoc(2);  % 纬度
        lon = geoLoc(1);  % 经度
        timezone = 8;     % UTC+8 不同地区时区不一样，中国地区统一采用+8
        % 循环处理每个时间点
        for i = 1:length(solarOpt.time)
            % 解析当前时间
            [year, month, day, hour, ~, ~] = datevec(solarOpt.time(i));
            
            % 计算太阳天顶角
            sza = calculate_sza(lat, lon, year, month, day, hour , timezone);
            
            % 自定义SBDART输入文件内容
            input_content = {
                '&INPUT'
                'iout  = 1,'
                sprintf('wlinf = %.3f,', p.Results.spectrumLimits(1)/1000)  %  起始步长
                sprintf('wlsup = %.3f,', p.Results.spectrumLimits(2)/1000)  %  终止步长
                sprintf('wlinc = %.3f,', 5/1000)  % 固定步长
                'idatm = 3,'    % 中纬度冬季大气
                'iaer  = 2,'    % 城市气溶胶
                'zbaer  = 1.5, 3.0,'
                'dbaer  = 0.3, 0.1, '
                'tbaer  = 0.8, '
                'isalb = 10,'    % 地表反照率
                'sc(1)  = 0.0,'  % 0%雪地+10%水体+20%沙地+70%植被
                'sc(2)  = 0.1,'
                'sc(3)  = 0.2,'
                'sc(4)  = 0.7,'
                'nstr  = 16,'    % 提升流数以提高辐射率精度、
                'tcloud = 20,1'  % 云光学厚度（层积云典型值：10-30）
                'zcloud = 1,-2'   % 云底高度（km，武汉层积云常见高度1-2km）
                sprintf('sza   = %.2f,', sza)
                'phi   = 0,'        % 方位角（0°为太阳方向）
                '/'
            };
            
            % 写入输入文件
            fid = fopen(fullfile(sbdartInputDir, 'INPUT'), 'w');
            fprintf(fid, '%s\n', input_content{:});
            fclose(fid);
            
            % 执行SBDART命令
            cygwin_command = sprintf('cd "/cygdrive/c/Users/dong/Desktop/SBDART-master" && %s', ...
                cygwinSbdartPath);
            [status, cmdout] = system(sprintf('"%s" -l -c "%s"', sbdartPath, cygwin_command));
            
            if status ~= 0
                error('SBDART执行失败（时间点 %d）: %s', i, cmdout);
            end
            
            % 解析输出数据
            outputLines = strsplit(cmdout, '\n');
            dataStart = 4;  % 跳过3行标题
            validDataLines = outputLines(dataStart:end-1);
            data = textscan(strjoin(validDataLines, '\n'), '%f %f %f %f %f %f %f %f');
            
            % 提取直射辐射（第6列）并匹配波长 W/m²/µm → W/m²/nm/sr^-1
            irradianceData = data{6}/1000/pi; % sky spectral irradiance (W m^-2 nm^-1 sr^-1)
            if length(irradianceData) == length(solarOpt.wavelength)
                solarOpt.irradiance(i,: ) = irradianceData;  
            else
                error('波长数据维度不匹配（时间点 %d）: 期望 %d，实际 %d',...
                      i, length(solarOpt.wavelength), length(irradianceData));
            end
        end
        
    case 'daytime'
        skyDiffusedDatabase = fullfile(LISAR_ENVS.RootDir, 'data', 'sky-diffused-irradiance.mat');
        data = load(skyDiffusedDatabase, 'wavelength', 'irradiance');
        isNotNan = (~ isnan(data.wavelength)) & (~ isnan(data.irradiance));
        data.wavelength = data.wavelength(isNotNan);
        data.irradiance = data.irradiance(isNotNan);
        solarOpt.wavelength = p.Results.spectrumLimits(1):p.Results.spectrumLimits(2);
        solarOpt.irradiance = repmat(interp1(data.wavelength, data.irradiance, solarOpt.wavelength, 'linear', 0) / (4 * pi), length(solarOpt.time), 1);

    case 'nighttime'
        solarOpt.wavelength = p.Results.spectrumLimits(1):p.Results.spectrumLimits(2);
        solarOpt.irradiance = zeros(length(solarOpt.time), length(solarOpt.wavelength));

    otherwise
        error('Unknown solar sky irradiance model: %s', p.Results.model);
end

end