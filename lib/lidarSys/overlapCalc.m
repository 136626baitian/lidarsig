function [Oz] = overlapCalc(height, theta_laser, theta_tel, diameter_tel, xOff, yOff, xMisAlign, yMisAlign, varargin)
% OVERLAPCALC Calculate overlap function of lidar.
%
% USAGE:
%    [Oz] = overlapCalc(height, theta_laser, theta_tel, xOff, yOff, xMisAlign, yMisAlign)
%
% INPUTS:
%    height: numeric
%        distance to the lidar system along line-of-sight. (m) 
%    theta_laser: numeric
%        divergence of laser beam. (rad)
%    theta_tel: numeric
%        field of fiew. (rad)
%    diameter_tel: numeric
%        diameter of the telescope. (m)
%    xOff: numeric
%        offset along x-axis. (m)
%    yOff: numeric
%        offset along y-axis. (m)
%    xMisAlign: numeric
%        misalignment along x-axis. (rad)
%    yMisAlign: numeric
%        misalignment along y-axis. (rad)
%
% OUTPUTS:
%    Oz: numeric
%        overlap function.
%
% HISTORY:
%    2022-06-15: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'height', @isnumeric);
addRequired(p, 'theta_laser', @isnumeric);
addRequired(p, 'theta_tel', @isnumeric);
addRequired(p, 'diameter_tel', @isnumeric);
addRequired(p, 'xOff', @isnumeric);
addRequired(p, 'yOff', @isnumeric);
addRequired(p, 'xMisAlign', @isnumeric);
addRequired(p, 'yMisAlign', @isnumeric);

parse(p, height, theta_laser, theta_tel, diameter_tel, xOff, yOff, xMisAlign, yMisAlign, varargin{:});

Oz = zeros(size(height));

for iH = 1:length(height)
    wr = theta_laser / 2 * height(iH);
    d = sqrt((xOff + xMisAlign * height(iH) - diameter_tel / 2).^2 + (yOff + yMisAlign * height(iH)).^2);
    func = @(r, phi)(2 ./ (pi * wr.^2) .* exp(-2 * (r.^2 + d.^2 - 2 .* r .* d .* cos(phi)) ./ wr.^2) .* r);
    Oz(iH) = integral2(func, 0, theta_tel / 2 * height(iH), 0, 2 * pi);
end

end