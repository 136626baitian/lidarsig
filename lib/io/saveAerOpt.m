function saveAerOpt(filepath, time, height, wavelength, extinction, backscatter, depolarization_ratio)
% SAVEAEROPT save profiles of aerosol/cloud optical property
%
% USAGE:
%    saveAerOpt(filepath, time, height, wavelength, extinction, backscatter, depolarization_ratio)
%
% INPUTS:
%    filepath: char
%        absolute path of output netcdf4 file.
%    time: array
%        matlab time stamp.
%    height: array
%        height array (m)
%    wavelength: array
%        incident laser wavelength. (nm)
%    extinction: matrix (time x height x wavelength)
%        extinction coefficient. (m-1)
%    backscatter: matrix (time x height x wavelength)
%        backscatter coefficient. (m-1 sr-1)
%    depolarization_ratio: matrix (time x height x wavelength)
%        depolarization ratio.
%
% EXAMPLE:
%
% HISTORY:
%    2024-08-08: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode, netcdf.getConstant('CLASSIC_MODEL'));
mode = bitor(mode, netcdf.getConstant('CLOBBER'));
ncID = netcdf.create(filepath, mode);

% define dimensions
dimID_height = netcdf.defDim(ncID, 'height', length(height));
dimID_time = netcdf.defDim(ncID, 'time', length(time));
dimID_wavelength = netcdf.defDim(ncID, 'wavelength', length(wavelength));

% define variables
varID_height = netcdf.defVar(ncID, 'height', 'NC_FLOAT', dimID_height);
varID_time = netcdf.defVar(ncID, 'time', 'NC_FLOAT', dimID_time);
varID_wavelength = netcdf.defVar(ncID, 'wavelength', 'NC_FLOAT', dimID_wavelength);
varID_extinction = netcdf.defVar(ncID, 'extinction', 'NC_FLOAT', [dimID_time, dimID_height, dimID_wavelength]);
varID_backscatter = netcdf.defVar(ncID, 'backscatter', 'NC_FLOAT', [dimID_time, dimID_height, dimID_wavelength]);
varID_depolarization = netcdf.defVar(ncID, 'depolarization_ratio', 'NC_FLOAT', [dimID_time, dimID_height, dimID_wavelength]);

% define filling values
netcdf.defVarFill(ncID, varID_extinction, false, -9999);
netcdf.defVarFill(ncID, varID_backscatter, false, -9999);
netcdf.defVarFill(ncID, varID_depolarization, false, -9999);

% define data compression
netcdf.defVarDeflate(ncID, varID_extinction, true, true, 5);
netcdf.defVarDeflate(ncID, varID_backscatter, true, true, 5);
netcdf.defVarDeflate(ncID, varID_depolarization, true, true, 5);

% leave define mode
netcdf.endDef(ncID);

% write data to NC file
netcdf.putVar(ncID, varID_height, height);
netcdf.putVar(ncID, varID_time, datenum_2_unix_timestamp(time));
netcdf.putVar(ncID, varID_wavelength, wavelength);
netcdf.putVar(ncID, varID_extinction, extinction);
netcdf.putVar(ncID, varID_backscatter, backscatter);
netcdf.putVar(ncID, varID_depolarization, depolarization_ratio);

% re-enter define mode
netcdf.reDef(ncID);

netcdf.putAtt(ncID, varID_height, 'unit', 'm');
netcdf.putAtt(ncID, varID_height, 'long_name', 'height above sea level');
netcdf.putAtt(ncID, varID_height, 'standard_name', 'height');

netcdf.putAtt(ncID, varID_time, 'unit', 'seconds since 1970-01-01 00:00:00 UTC');
netcdf.putAtt(ncID, varID_time, 'long_name', 'Time UTC');
netcdf.putAtt(ncID, varID_time, 'standard_name', 'time');
netcdf.putAtt(ncID, varID_time, 'axis', 'T');
netcdf.putAtt(ncID, varID_time, 'calendar', 'julian');

netcdf.putAtt(ncID, varID_wavelength, 'unit', 'nm');
netcdf.putAtt(ncID, varID_wavelength, 'long_name', 'wavelength');
netcdf.putAtt(ncID, varID_wavelength, 'standard_name', 'wavelength');

netcdf.putAtt(ncID, varID_extinction, 'unit', 'm^-1');
netcdf.putAtt(ncID, varID_extinction, 'long_name', 'extinction coefficient');
netcdf.putAtt(ncID, varID_extinction, 'standard_name', 'alpha');

netcdf.putAtt(ncID, varID_backscatter, 'unit', 'm^-1 sr^-1');
netcdf.putAtt(ncID, varID_backscatter, 'long_name', 'backscatter coefficient');
netcdf.putAtt(ncID, varID_backscatter, 'standard_name', 'beta');

netcdf.putAtt(ncID, varID_depolarization, 'unit', '');
netcdf.putAtt(ncID, varID_depolarization, 'long_name', 'depolarization_ratio');
netcdf.putAtt(ncID, varID_depolarization, 'standard_name', 'delta');

% write global attributes
varID_global = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncID, varID_global, 'License', 'CF-1.0');

% close NC file
netcdf.close(ncID);

end