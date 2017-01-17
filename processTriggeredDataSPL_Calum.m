function varargout = processTriggeredDataSPL(varargin)
%------------------------------------------------------------------------
% dbvals = processTriggeredDataSPL(	
%										'base', <base file name for .daq files>,
% 										'path',	<path to .daq files>,
% 										'sense', <microphone sensitivity>, 
% 										'gain', <microphone gain>, 
% 										'window', <RMS window size>,
% 										'lpfc', <lowpass filter freq>,
% 										'lporder', <lowpass filter order>,
% 										'hpfc', <highpass filter freq>,
% 										'hporder', <highpass filter order> )
%------------------------------------------------------------------------
% TytoLogy:NICal program
%------------------------------------------------------------------------
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
% 
%------------------------------------------------------------------------
% See also: NICal 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 20 January, 2016 from processTriggeredDataRMSMJR.m (SJS)
%
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Variables and Constants Declarations
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% window size (in milliseconds) for computing rms (and dB) values
%	use smaller values for greater resolution, larger for coarse resolution
rms_windowsize_ms = 10;

% highpass cutoff frequency (Hz)
hpFc = 1000;
% filter order
hporder = 3;
% lowpass cutoff Frequency (Hz)
lpFc = 200000;
lporder = 3;

% Microphone volts/Pa and gain
MicSense = 1;
MicGain = 1;
VtoPa = (1./MicSense)*(1./MicGain);	

% Decimation factor - plotted data will be 1 / DeciFactor shorter
% and sampling rate will be Fs / DeciFactor
DeciFactor = 1;

% set basepath and basename to empty and
% use default calibration mode ('tones')
basepath = '';
basename = '';

% Get input from user in case signal is offset from trigger
% winstart = input('Enter any delay to start of signal in ms: ');
% if isempty(winstart)
% 	winstart=0;
% end

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
			
		case 'sense'
			if ~isnumeric(varargin{index+1})
				error('%s: mic sense value must be a number', mfilename);
			elseif ~isscalar(varargin{index+1})
				error('%s: mic sense value must be a single value', mfilename);
			elseif varargin{index+1} <= 0
				error('%s: mic sense value must be > 0', mfilename);
			else
				MicSense = varargin{index+1};
			end
			index = index + 2;
			
		case 'gain'
			if ~isnumeric(varargin{index+1})
				error('%s: gain value must be a number', mfilename);
			elseif ~isscalar(varargin{index+1})
				error('%s: gain value must be a single value', mfilename);
			elseif varargin{index+1} <= 0
				error('%s: gain value must be > 0', mfilename);
			else
				MicGain = varargin{index+1};
			end
			index = index + 2;

		case 'window'
			if ~isnumeric(varargin{index+1})
				error('%s: window value must be a number', mfilename);
			elseif ~isscalar(varargin{index+1})
				error('%s: window value must be a single value', mfilename);
			elseif varargin{index+1} <= 0
				error('%s: window value must be > 0', mfilename);
			else
				rms_windowsize_ms = varargin{index+1};
			end
			index = index + 2;
	
		case 'hpfc'
			if ~isnumeric(varargin{index+1})
				error('%s: highpass Fc value must be a number', mfilename);
			elseif ~isscalar(varargin{index+1})
				error('%s: highpass Fc value must be a single value', mfilename);
			elseif varargin{index+1} <= 0
				error('%s: highpass Fc value must be > 0', mfilename);
			else
				hpFc = varargin{index+1};
			end
			index = index + 2;
			
		case 'hporder'
			if ~isnumeric(varargin{index+1})
				error('%s: highpass order value must be a number', mfilename);
			elseif ~isscalar(varargin{index+1})
				error('%s: highpass order value must be a single value', mfilename);
			elseif varargin{index+1} <= 0
				error('%s: highpass order value must be > 0', mfilename);
			else
				hporder = varargin{index+1};
			end
			index = index + 2;			
			
		case 'lpfc'
			if ~isnumeric(varargin{index+1})
				error('%s: lowpass Fc value must be a number', mfilename);
			elseif ~isscalar(varargin{index+1})
				error('%s: lowpass Fc value must be a single value', mfilename);
			elseif varargin{index+1} <= 0
				error('%s: lowpass Fc value must be > 0', mfilename);
			else
				lpFc = varargin{index+1};
			end
			index = index + 2;
			
		case 'lporder'
			if ~isnumeric(varargin{index+1})
				error('%s: lowpass order value must be a number', mfilename);
			elseif ~isscalar(varargin{index+1})
				error('%s: lowpass order value must be a single value', mfilename);
			elseif varargin{index+1} <= 0
				error('%s: lowpass order value must be > 0', mfilename);
			else
				lporder = varargin{index+1};
			end
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
	% if not, ask user for file
	
	%----------------------------------
	% output data path and file
	%----------------------------------
	% DefaultPath = 'C:\Users\Calibrate\Matlab\Calibration\Data';
	DefaultPath = pwd;
	
	%----------------------------------
	% open panel to get .daq file name
	%----------------------------------
	[tmpfile, basepath] = uigetfile(...
			 {'*.daq', 'DAQ output files (*.daq)'}, ...
			  'Pick a DAQ file', DefaultPath);
	% check if user hit cancel (tmpfile, or basepath == 0)
	if isequal(tmpfile, 0) || isequal(basepath, 0) || isempty(tmpfile)
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

%-----------------------------------------------
% find the daq files that match the basename
%-----------------------------------------------
[daqFiles, daqNumbers] = find_daq_files(basename, basepath);
nDaqFiles = length(daqFiles);
% error if # daq files is 0
if ~nDaqFiles
	error('%s: no daq files found (%s, %s)', mfilename, basename, basepath);
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
rmsvals = cell(nDaqFiles, 1);
dbvals = cell(nDaqFiles, 1);

for n = 1:nDaqFiles

	%--------------------------------
	% build filename
	%--------------------------------
	daqfile = fullfile(basepath, daqFiles{n});
	
	%--------------------------------
	% Read Data
	%--------------------------------
	[tmpdata, t, tabs, events, info] = daqread(daqfile);
	[npts, nchan] = size(tmpdata);
	fprintf('Read %d points from file %s\n', npts, daqfile);
	if isempty(tmpdata)
		% file is empty
		fprintf('File %s is empty, terminating data loading\n\n', daqfile);
		break
	end	
	% ASSUME (!!) that stimulus data are collected on channel 1 (NI 0)
	% and mic data are on channel 2 (NI AI0)
	stimdata = tmpdata(:, 1);
	micdata = tmpdata(:, 1);
	clear tmpdata
	% get sample rate from the DAQ info object
	Fs = info.ObjInfo.SampleRate;

	%------------------------------------------------------------------------
	% get a highpass and lowpass filter for processing the  data
	%------------------------------------------------------------------------
	% Nyquist frequency
	fnyq = Fs/2;
	% filter coefficients
	[HPb, HPa] = butter(hporder, hpFc/fnyq, 'high');
	[LPb, LPa] = butter(lporder, lpFc/fnyq, 'low');

	%------------------------------------------------------------------------
	% now filter data
	%------------------------------------------------------------------------
	% window and filter the data
	micdata = sin2array(micdata', 2, Fs);
	micdata = filter(HPb, HPa, micdata);
	micdata = filter(LPb, LPa, micdata);
	
	%--------------------------------
	% process data 
	%--------------------------------
	[rmsvals{n}, rmw_windows{n}] =  ...
									processWindows(micdata, rms_windowsize_ms, Fs);
	% convert to dB SPL
	dbvals{n} = dbspl(VtoPa * rmsvals{n});

	% decimate data for plotting
	micdata_reduced = decimate(micdata, DeciFactor);
	Fs_reduced = Fs / 10;
	% build time vectors for plotting
	t1 = ((1:length(micdata_reduced)) - 1) / Fs_reduced;
	t2 = rms_windowsize_ms * 0.001 * (0:size(dbvals{n},1)-1);
% 			t2 = rms_windowsize_ms * 0.001 * (0:rmsIndex);
	% plot!

	figure(n)
	subplot(211)
	plot(t1, micdata_reduced);
	title(sprintf('File: %s', daqFiles{n}));
	grid
	ylabel('Volts');
	subplot(212)
	plot(t2, dbvals{n}, 'Marker', '.', 'Color', 'r');
	ylabel('dB SPL')
	xlabel('Time (seconds)');
	ylim([0 100]);
	grid
	drawnow
end

% assign output vars
out.dbvals = dbvals;
out.rmsvals = rmsvals;
out.rms_windowsize_ms = rms_windowsize_ms;

out.Fs = Fs;
out.files = daqFiles;
out.path = basepath;
out.HP = struct('fc', hpFc, 'order', hporder, 'a', HPa, 'b', HPb);
out.LP = struct('fc', lpFc, 'order', lporder, 'a', LPa, 'b', LPb);
out.DeciFactor = DeciFactor;
out.VtoPa = VtoPa;
varargout{1} = out;

end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% find_daq_files function
%------------------------------------------------------------------------
%------------------------------------------------------------------------
function varargout = find_daq_files(basename, varargin)
	if nargin == 2
		basepath = varargin{1};
	else
		basepath = pwd;
	end

	% build search string for dir command from the basename
	searchStr = sprintf('%s-*.daq', basename);
	% get struct array from dir command
	tmpstruct = dir(fullfile(basepath, searchStr));
	% convert to cell array
	tmpcell = struct2cell(tmpstruct);
	% keep the names
	daqFileNames = tmpcell(1, :);
	% clear unused variables
	clear tmpstruct tmpcell;

	% get # of found daqFileNames
	nDaqFiles = length(daqFileNames);
	if nDaqFiles == 0
		varargout{1} = '';
		varargout{2} = [];
		return
	end

	% allocate vector for daq file numbers
	daqFileNumbers = zeros(nDaqFiles, 1);
	% loop through the daq file names
	for n = 1:nDaqFiles
		% get the filename without the extension
		[~, tmpname, ~] = fileparts(daqFileNames{n});
		% find the dashes in the filename
		dashindices = find(tmpname == '-');
		% take the characters after the final dash and convert to number,
		% then store this number in daqFileNumbers
		daqFileNumbers(n) = str2num(tmpname( (dashindices(end)+1):end)); %#ok<ST2NM>
	end

	% sort the numbers and reorder the daqFileNames
	[~, sorted_indices] = sort(daqFileNumbers);
	varargout{2} = daqFileNumbers(sorted_indices);
	varargout{1} = daqFileNames(sorted_indices);
end




%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% processWindows function
%------------------------------------------------------------------------
%------------------------------------------------------------------------
function [rmsvals, rmswindows] = processWindows(data, windowms, fs)
	% convert windowsize from msec into samples
	windowpts = ms2samples(windowms, fs);
	
	% calculate rms_window indices into data
	rmswindows = 1:windowpts:length(data);
	Nwindows = length(rmswindows);
	
	% compute rms values for all windows of data
	rmsvals = zeros(Nwindows, 1);
	rmsIndex = 0;
	for w = 2:Nwindows
		rmsIndex = rmsIndex + 1;
		rmsvals(rmsIndex) = rms(data(rmswindows(w-1):rmswindows(w)));
	end
end
