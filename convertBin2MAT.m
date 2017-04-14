function varargout = convertBin2MAT(varargin)
%------------------------------------------------------------------------
% D = convertTrig2MAT(	'inputfile', <input .bin file>,
% 								'outputfile', <output path/file', 
%)
%------------------------------------------------------------------------
% TytoLogy:NICal program
%------------------------------------------------------------------------
% FOR SESSION DATA ONLY!!!!!!!
% 
% If inputfile is provided, program will read raw data from that file
% Otherwise, user will be prompted for input file
% 
% If outputfile is provided, output will be stored in specified file.
% Otherwise, output will be stored in file <basename>_convert.mat
% 
% E.g., 
% 
% convertBin2MAT(	'inputfile', 'data.bin', ...
% 						'outputfile', 'D:\calibrationuser\data.mat' )
%
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
inputfile = ''; %#ok<NASGU>
outputfile = ''; %#ok<NASGU>
matfile = ''; %#ok<NASGU>

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
		case 'inputfile'
			% set inputfile
			inputfile = varargin{index+1}; %#ok<NASGU>
			% increment index by 2 places
			index = index + 2;
			
		case 'outputfile'
			outputfile = varargin{index+1}; %#ok<NASGU>
			% increment index by 2 places
			index = index + 2;
		
		otherwise
			error('%s: invalid option %s', mfilename, varargin{index});
	end	
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Variables and Constants Declarations
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set basepath and basename to empty and
% use default calibration mode ('tones')
inputfile = 'optotest_1.bin';
outputfile = 'optotest_1_converted.mat';

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% find data files
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%----------------------------------
% check if input filename was provided
%----------------------------------
if isempty(inputfile)
	% if so, ask user for file
	
	%----------------------------------
	% output data path and file
	%----------------------------------
	DefaultPath = pwd;
	
	%----------------------------------
	% open panel to get .daq file name
	%----------------------------------
	[inputfile, basepath] = uigetfile(...
			 {'*.bin', 'BIN output files (*.bin)'}, ...
			  'Pick a .BIN file', DefaultPath);
	% check if user hit cancel (tmpfile, or basepath == 0)
	if isequal(inputfile, 0)
		disp('Cancelled...')
		varargout{1} = [];
		return
	else
		inputfile = fullfile(basepath, inputfile);
	end
	
else
	% inputfile was given, see if it exists 
	if ~exist(inputfile, 'file')
		error('%s: inputfile %s not found', mfilename, inputfile)
	end
end

%----------------------------------
% break down file into parts
%----------------------------------
[basepath, basename, ~] = fileparts(inputfile);

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
	if isequal(outputfile, 0)
		disp('Cancelled...')
		return
	else
		% combine filename and path
		outputfile = fullfile(outputpath, outputfile);
	end
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% check for matfile
%------------------------------------------------------------------------
%------------------------------------------------------------------------
matfile = fullfile(basepath, [basename '.mat']);
if ~exist(matfile, 'file')
	warning('%s: cannot locate NICal information file %s', mfilename, matfile);
	cal = [];
else
	load(matfile, '-MAT');
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% READ DATA
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%-----------------------------------------------
% open file
%-----------------------------------------------
fp = fopen(inputfile, 'r');

%-----------------------------------------------
% read cal struct
%-----------------------------------------------
tmp = readStruct(fp);
% if cal is empty (matfile not found), get cal information from file
if isempty(cal)
	cal = tmp;
	% convert character fields
	cal.mic_fr_file = char(cal.mic_fr_file);
	cal.calfile = char(cal.calfile);
	cal.TriggerSettings.TriggerType = char(cal.TriggerSettings.TriggerType);
	cal.TriggerSettings.TriggerSource = char(cal.TriggerSettings.TriggerSource);
	cal.TriggerSettings.TriggerCondition = ...
										char(cal.TriggerSettings.TriggerCondition);
end
clear tmp

%-----------------------------------------------
% read sweeps
%-----------------------------------------------
% allocate temporary storage
if any(cal.Side == [1 2])
	tmp = cell(1000, 1);
else
	tmp = cell(1000, 2);	
end

fFlag = 0;
index = 1;
while(~feof(fp) && ~fFlag)
	% read L data
	if feof(fp)
		fFlag = 1;
	else
		if any(cal.Side == [1 3])
			tmp{index, 1} = readVector(fp);
		end
	end
	
	% read R data
	if feof(fp)
		fFlag = 1;
	else
		if any(cal.Side == [2 3])
			tmp{index, 2} = readVector(fp);
		end
	end
	
	% check for eof or error
	[fmsg, ferr] = ferror(fp);
	if feof(fp) || ferr
		fFlag = 1;
		fprintf('%s: end of file... message: <%s>\n', mfilename, fmsg);
	else
		index = index + 1;
	end
end

%-----------------------------------------------
% close file
%-----------------------------------------------
fclose(fp);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% assign to D struct
%------------------------------------------------------------------------
%------------------------------------------------------------------------
D.cal = cal;
% check that last sweep is not empty
if (isempty(tmp{index, 1}) && isempty(tmp{index, 2}))
	% if so, move back 1
	index = index - 1;
	if index < 1
		error('%s: odd, or empty data file %s', mfilename, inputfile);
	end
end
D.data = tmp(1:index, :);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% save to output file
%------------------------------------------------------------------------
%------------------------------------------------------------------------
save(outputfile, 'D', '-MAT');

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% assign output vars
%------------------------------------------------------------------------
%------------------------------------------------------------------------
varargout{1} = D;
