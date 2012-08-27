function varargout = processNICal(varargin)
%------------------------------------------------------------------------
% dbvals = processDAQdB(basename, basepath, calmode)
%------------------------------------------------------------------------
% 
% 	If basename and basepath are provided, program will use those values 
% 	to search for .daq files.
% 	
% 	Otherwise, a GUI dialog will pop up to ask user for file
% 	
% 	calmode = 'tones'  (default) will look for tone magnitudes in
% 					each .daq file and compute a freq-dB SPL curve
% 					
% 	calmode = 'rms'	will simply compute the overall db SPL level
% 							for each file
%	calmode = 'window' will compute db SPL levels for each file 
% 							broken up into 100 msec windows
% 								
%------------------------------------------------------------------------
% See also: DAQ Toolbox 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 29 February, 2012 (SJS)
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
rms_windowsize_ms = 100;

% highpass cutoff frequency (Hz)
fcutoff = 90;
% filter order
forder = 3;

% Decimation factor - plotted data will be 1 / DeciFactor shorter
% and sampling rate will be Fs / DeciFactor
DeciFactor = 10;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% parse inputs
%------------------------------------------------------------------------
%------------------------------------------------------------------------

if nargin == 0
	basepath = '';
	basename = '';
	calmode = 'tones';
elseif nargin == 1
	[~, basename] = fileparts(varargin{1});
	basepath = pwd;
	calmode = 'tones';
elseif nargin == 2
	[~, basename] = fileparts(varargin{1});
	basepath = varargin{2};
	calmode = 'tones';
elseif nargin == 3
	[~, basename] = fileparts(varargin{1});
	basepath = varargin{2};
	calmode = varargin{3};
	if ~any(strcmpi(calmode, {'tones', 'rms', 'window'}))
		error('%s: invalid calmode %s', mfilename, calmode);
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
	micdata = tmpdata(:, 2);
	clear tmpdata
	
	% get sample rate
	Fs = info.ObjInfo.SampleRate;

	%------------------------------------------------------------------------
	% get a highpass filter for processing the  data
	%------------------------------------------------------------------------
	% Nyquist frequency
	fnyq = Fs/2;
	% filter coefficients
	[fcoeffb, fcoeffa] = butter(forder, fcutoff/fnyq, 'high');

	%------------------------------------------------------------------------
	% now process data
	%------------------------------------------------------------------------

	% window and filter the data
	micdata = sin2array(micdata', 1, Fs);
	micdata = filter(fcoeffb, fcoeffa, micdata);
	
	switch lower(calmode)
		case 'window'
			% convert rms_windowsize from msec into samples
			rms_windowsize = ms2samples(rms_windowsize_ms, Fs);
			% calculate rms_window indices into data
			rms_windows = 1:rms_windowsize:length(micdata);
			Nwindows = length(rms_windows);
			% compute rms values for all windows of data
			rms_vals = zeros(Nwindows, 1);
			rmsIndex = 0;
			for w = 2:Nwindows
				rmsIndex = rmsIndex + 1;
				rms_vals(rmsIndex) = rms(micdata(rms_windows(w-1):rms_windows(w)));
			end
			% convert to dB SPL
			dbvals{n} = dbspl(VtoPa * rms_vals);
			% decimate data for plotting
			micdata_reduced = decimate(micdata, DeciFactor);
			Fs_reduced = Fs / 10;
			% build time vectors for plotting
			t1 = ((1:length(micdata_reduced)) - 1) / Fs_reduced;
			t2 = rms_windowsize_ms * 0.001 * (0:rmsIndex);
			% plot!

			subplot(211)
			plot(t1, micdata_reduced);
			grid
			ylabel('Volts');

			subplot(212)
			plot(t2, dbvals{n}, 'Marker', '.', 'Color', 'r');
			ylabel('dB SPL')
			xlabel('Time (seconds)');
			grid

		case 'tones'
			[tmpfreqs, tmpmags, fmax, magmax] = daqdbfft(micdata, Fs, length(micdata));
			calfreqs(n) = fmax;
			[mags(n), phis(n)] = fitsinvec(micdata, 1, Fs, fmax);
	end
	
end

switch lower(calmode)
	case 'tones'
		ndatums = length(calfreqs);
		out.freqs = calfreqs;
		out.mags = mags;
		out.phis = phis;
		out.dbvals = dbspl(VtoPa * rmssin * out.mags);

		figure
		subplot(211)
		plot(dbspl(VtoPa * rmssin * out.mags), '.-')
		set(gca, 'XTick', 1:ndatums);
		set(gca, 'XTickLabel', '');
		xlim([0 ndatums+1])
		grid
		title(basename);
		ylabel('dB SPL')

		subplot(212)
		plot(unwrap(out.phis), '.-')
		set(gca, 'XTick', 1:ndatums);
		for n = 1:ndatums
			labels{n} = sprintf('%.1f', 0.001*calfreqs(n));
		end
		set(gca, 'XTickLabel', labels);
		xlim([0 ndatums+1])
		grid
		ylabel('phase (rad)');
		xlabel('frequency (kHz)');

	case 'window'
		out.dbvals = dbvals;
		out.rms_windowsize_ms = rms_windowsize_ms;

end

out.Fs = Fs;
out.files = daqFiles;
out.path = basepath;
out.fcutoff = fcutoff;
out.forder = forder;
out.DeciFactor = DeciFactor;
out.fcoeff.a = fcoeffa;
out.fcoeff.b = fcoeffb;

varargout{1} = out;