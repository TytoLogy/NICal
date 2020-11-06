function varargout = convertTrig2MAT(varargin)
%------------------------------------------------------------------------
% dbvals = convertTrig2MAT(	'base', <base file name for .daq files>,
% 										'path',	<path to .daq files>, 
% 										'outputfile', <output path/file')
%------------------------------------------------------------------------
% TytoLogy:NICal program
%------------------------------------------------------------------------
% NOT FOR SESSION DATA!!!!!!!
% 
% If base filename and base path are provided, program will use those values 
% to search for .daq files.
% 
% The base filename should be the root filename of the automatically
% triggered .daq files from NICal.  For example, if the files collected by
% NICal are:
% 	NICalfile_test-1.daq
% 	NICalfile_test-2.daq
% 	.
% 	.
% 	.
% 	NICalfile_test-n.daq
% and the files are in the directory D:\calibrationuser\data, base and path 
% would be declared as follows:
% 
% processTriggeredData(	'base', 'NICalfile_test', ...
% 								'path', 'D:\calibrationuser\data' )
% 
% If basename or path is not provided, a GUI dialog will pop up
% to ask user for file
%
% If outputfile is provided, output will be stored in specified file.
% Otherwise, output will be stored in file <basename>_convert.mat
% 
%------------------------------------------------------------------------
% See also: NICal 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 14 April, 2017 from processTriggeredData (SJS)
%
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Variables and Constants Declarations
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set basepath and basename to empty and
% use default calibration mode ('tones')
basepath = '';
basename = '';
outputfile = '';

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% parse inputs
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% loop through arguments
index = 1;
while index <= nargin
	% act according to tag
	switch lower(varargin{index})
		case 'base'
			% set basename, strip off any path information
			[~, basename] = fileparts(varargin{index+1});
			% increment index by 2 places
			index = index + 2;

		case 'path'
			basepath = varargin{index+1};
			% increment index by 2 places
			index = index + 2;
			
		case 'outputfile'
			outputfile = varargin{index+1};
			% increment index by 2 places
			index = index + 2;
		
		otherwise
			error('%s: invalid option %s', mfilename, varargin{index});
	end	
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% find data files
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%----------------------------------
% check if filename was provided
%----------------------------------
if isempty(basename)
	% if so, ask user for file
	
	%----------------------------------
	% output data path and file
	%----------------------------------
	DefaultPath = pwd;
	
	%----------------------------------
	% open panel to get .daq file name
	%----------------------------------
	[tmpfile, basepath] = uigetfile(...
			 {'*.daq', 'DAQ output files (*.daq)'}, ...
			  'Pick a DAQ file', DefaultPath);
	% check if user hit cancel (tmpfile, or basepath == 0)
	if isequal(tmpfile, 0) || isequal(basepath, 0)
		disp('Cancelled...')
		return
	end
	
	%----------------------------------
	% break down file into parts
	%----------------------------------
	[~, extfile, ~] = fileparts(tmpfile);
	% last 2 chars of file should be a '-' and a number; remove this to get
	% a valid basename for file
	basename = extfile(1:(length(extfile) - 2));
	
else
	% basename was given, see if basepath was also provided
	if ~exist('basepath', 'var') || isempty(basepath)
		% if not, use current dir
		basepath = pwd;
	end
end

%----------------------------------
% check if outputfile was provided
%----------------------------------
if isempty(outputfile)
	% if so, ask user for file
	%----------------------------------
	% output mat path and file
	%----------------------------------
	DefaultPath = pwd;
	DefaultFile = [basename '_convert.mat'];
	
	%----------------------------------
	% open panel to get .mat file name
	%----------------------------------
	[outputfile, outputpath] = uiputfile(...
			 {'*.mat', 'MAT-files (*.mat)'}, ...
			  'Save As', fullfile(DefaultPath, DefaultFile) );
	% check if user hit cancel (tmpfile, or basepath == 0)
	if isequal(tmpfile, 0)
		disp('Cancelled...')
		return
	else
		outputfile = fullfile(outputpath, outputfile);
	end
end

%-----------------------------------------------
% find the daq files that match the basename
%-----------------------------------------------
[daqFiles, ~] = find_daq_files(basename, basepath);
nDaqFiles = length(daqFiles);
% error if # daq files is 0
if ~nDaqFiles
	error('%s: no daq files found (%s, %s)', mfilename, basename, basepath);
else
	% preallocate output struct
	D = repmat( struct(	'filename', '', ...
								'data', [], ...
								'info', [], ...
								'events', [] ), ...
					nDaqFiles, 1);
end

%-----------------------------------------------
% load .mat file for this calibration session
%-----------------------------------------------
load(fullfile(basepath, [basename '.mat']), '-MAT');

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% PROCESS DATA
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%--------------------------------
% loop through daqfiles
%--------------------------------
for n = 1:nDaqFiles

	%--------------------------------
	% build filename
	%--------------------------------
	daqfile = fullfile(basepath, daqFiles{n});
	
	%--------------------------------
	% Read Data
	%--------------------------------
	fprintf('Reading file %s...\n', daqfile);
	[tmpdata, t, tabs, events, info] = daqread(daqfile);
	[npts, nchan] = size(tmpdata);
	fprintf('Read %d points from file %s\n', npts, daqfile);
	if isempty(tmpdata)
		% file is empty
		fprintf('File %s is empty, terminating data loading\n\n', daqfile);
		break
	else
		% store data and information
		D(n).filename = daqfile;
		D(n).data = tmpdata;
		D(n).info = info;
		D(n).events = events;
	end
	clear tmpdata
end

% save to output file
save(outputfile, 'D', '-MAT');

% assign output vars
varargout{1} = D;

end

