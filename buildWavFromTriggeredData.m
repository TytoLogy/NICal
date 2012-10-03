function varargout = buildCalFromTriggeredData(varargin)
%------------------------------------------------------------------------
% cal = buildCalFromTriggeredData(basename, basepath, channel)
%------------------------------------------------------------------------
% TytoLogy:NICal program
%------------------------------------------------------------------------
% 
% 	If basename and basepath are provided, program will use those values 
% 	to search for .daq files.
% 	
% 	Otherwise, a GUI dialog will pop up to ask user for file
% 	
% 								
%------------------------------------------------------------------------
% See also: NICal, processTriggeredData, buildCalFromTriggeredData
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 3 October, 2012 (SJS) from buildCalFromTriggeredData
%
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Variables and Constants Declarations
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% highpass cutoff frequency (Hz)
fcutoff = 90;
% filter order
forder = 3;

% Decimation factor - plotted data will be 1 / DeciFactor shorter
% and sampling rate will be Fs / DeciFactor
DeciFactor = 10;
% AUTOFREQ set to 1 will use automagic tone freq detection
AUTOFREQ = 0;
% other variables that need declaration here
basepath = '';
basename = '';
channel = [];

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% parse inputs
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if nargin == 1
	[basepath, basename] = fileparts(varargin{1});
% 	basename = [basename '.daq'];
elseif nargin == 2
	[~, basename] = fileparts(varargin{1});
	basepath = varargin{2};
elseif nargin == 3
	[~, basename] = fileparts(varargin{1});
	basepath = varargin{2};
	channel = varargin{3};
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% check channels
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if isempty(channel)
	channel = query_userarray('Enter channel to convert', [], [1 2]);
end
	
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% check if filename was provided
%------------------------------------------------------------------------
%------------------------------------------------------------------------
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
	if isequal(tmpfile, 0) || isequal(basepath, 0)
		disp('Cancelled...')
		return
	end
	%----------------------------------
	% break down file into parts
	%----------------------------------
	[~, extfile, fext] = fileparts(tmpfile);
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

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% find the daq files that match the basename
%------------------------------------------------------------------------
%------------------------------------------------------------------------
[daqFiles, daqNumbers] = find_daq_files(basename, basepath);
nDaqFiles = length(daqFiles);
% error if # daq files is 0
if ~nDaqFiles
	error('%s: no daq files found (%s, %s)', mfilename, basename, basepath);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% load .mat file for this calibration session
%------------------------------------------------------------------------
%------------------------------------------------------------------------
load(fullfile(basepath, [basename '.mat']), '-MAT');

%--------------------------------------------------------
%--------------------------------------------------------
% PROCESS DATA
%--------------------------------------------------------
%--------------------------------------------------------
%--------------------------------
% loop through daqfiles
%--------------------------------
for n = 1:nDaqFiles

	%--------------------------------
	% build filename
	%--------------------------------
	daqfile = fullfile(basepath, daqFiles{n});
	[ptmp, btmp, extmp] = fileparts(daqfile);
	clear ptmp extmp
	
	%--------------------------------
	% Read Data
	%--------------------------------
	[tmpdata, t, tabs, events, info] = daqread(daqfile);
	[npts, nchan] = size(tmpdata);
	fprintf('Read %d points from file %s\n', npts, daqfile);
	if isempty(tmpdata)
		% file is empty
		fprintf('File %s is empty, terminating data loading\n\n', daqfile);
		nDaqFiles = nDaqFiles - 1;
		break
	end

	% get sample rate from the DAQ info object
	Fs = info.ObjInfo.SampleRate;

	%------------------------------------------------------------------------
	% write to .wav files
	%------------------------------------------------------------------------
	switch length(channel)
		case 1
			wfile = fullfile(basepath, sprintf('%s_ch%d.wav', btmp, channel));			
			wavwrite(0.99*normalize(tmpdata(:, channel)), Fs, 16, wfile);
		case 2
			for c = 1:length(channel)
				wfile = fullfile(basepath, sprintf('%s_ch%d.wav', btmp, channel(c)));			
				wavwrite(0.99*normalize(tmpdata(:, channel(c))), Fs, 16, wfile);
			end
		otherwise
			error('%s: channel is empty!', mfilename);
	end

end
