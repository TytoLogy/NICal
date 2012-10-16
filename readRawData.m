%--------------------------------------------------------------------------
% readRawData.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 16 October, 2012
% 
% Revisions:
%--------------------------------------------------------------------------

basepath = '/Users/sshanbhag/Work/Data/Audio/Calibration/MouseRig/10Oct2012';
basename = 'TDT7283_10Oct2012_45k-110k_4Vpk_FT6antialias'

datfile = fullfile(basepath, [basename, '.dat']);
calfile = fullfile(basepath, [basename, '.cal']);


%-----------------------------------------------------------------------
% check output  file - if it exists, check with user
%-----------------------------------------------------------------------
if exist(calfile, 'file')
	load(calfile, '-MAT', 'caldata')
else
	warning('NICAL:FileNotFound', '%s: could not find cal file', mfilename);
	fprintf('\t\tcal file: %s\n\n', calfile);

	[calfile, basepath] = uigetfile({'*.cal', '*_cal.mat'},'Open cal file', ...
													basepath);
	
	if isequal(calfile, 0) || isequal(basepath, 0)
		return
	else
		[~, basename, ~] = fileparts(calfile);
		calfile = fullfile(basepath, calfile);
		datfile = fullfile(basepath, [basename, '.dat']);
		load(calfile, '-MAT', 'caldata');
	end
end


