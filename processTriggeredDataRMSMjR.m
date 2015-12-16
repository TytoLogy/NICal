function varargout = processTriggeredDataRMSMjR(varargin)

clear all; close all;

%------------------------------------------------------------------------
% dbvals = processTriggeredData(	'base', <base file name for .daq files>,
% 											'path',	<path to .daq files>,
% 											'mode', <analysis mode>,
% 											'freqwidth' <automagic frequency detect width> )
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
% mode will select how the data are processed:
% 	mode = 'tones'  (default) will look for tone magnitudes in
% 					each .daq file and compute a freq-dB SPL curve
% 
% 	mode = 'rms'	will simply compute the overall db SPL level
% 							for each file
% 	mode = 'window' will compute db SPL levels for each file 
% 							broken up into 100 msec windows
% 
% freqwidth sets the width to search for peak frequency magnitude
% 		default is 21 Hz
%------------------------------------------------------------------------
% See also: NICal 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 29 February, 2012 (SJS)
%
% Revisions:
%	29 Aug 2012 (SJS): updated documentation
% 	22 Oct 2012 (SJS):
% 	 -	added freqdetectwidth input var
% 	 - reworking tone calibration
%	11 Dec 2012 (SJS):
%	 -	added section for rms calibration mode to output variable assignment
%	 - changed input parsing to pure varargin mode (no assumption on order 
%		of input args -> downside is that all inputs must have tags!
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

% tolerance (in Hz) for finding peak frequency in autodetect mode
FreqDetectWidth = 21;

% set basepath and basename to empty and
% use default calibration mode ('tones')
basepath = '';
basename = '';
calmode = 'rms';

% Get input from user in case signal is offset from trigger
winstart = input('Enter any delay to start of signal in ms: ');
if isempty(winstart)
	winstart=0;
end

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
			
		case 'mode'
			switch(lower(varargin{index+1}))
				case {'tones', 'window', 'rms'}
					calmode = lower(varargin{index+1});
				otherwise
					error('%s: invalid mode value %s', mfilename, varargin{index+1});
			end
			index = index + 2;

		case 'freqwidth'
			FreqDetectWidth = varargin{index+1};
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

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% should consider loading/asking user for 
% tone frequencies if tone is specified
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if strcmpi(calmode, 'tones')
	qVal = query_user('Auto-detect tone frequencies', 1);
	if qVal == 0
		% check if user has text list of frequencies
		fVal = query_user('Read list of frequencies from .txt file', 1);
		if fVal == 1
			% open panel to get .txt file name
			[freqfile, freqpath] = uigetfile(...
					 {'*.txt', 'txt frequency list files (*.txt)'}, ...
					  'Pick a .txt file');
			% check if user hit cancel (tmpfile, or basepath == 0)
			if isequal(freqfile, 0) || isequal(freqpath, 0)
				disp('Cancelled file load...')
				fVal = 0;
			else
				calfreqs = load(fullfile(freqpath, freqfile));
			end
		end
		% don't use else-if in order to allow fall-through from previous if
		if fVal == 0
			AUTOFREQ = 0;
			fstr = '';
			fstr = query_uservalue('Enter frequencies, separated by spaces', '');
			calfreqs = str2num(fstr);
			clear fstr;			
		end
	else
		AUTOFREQ = 1;
		calfreqs = zeros(1, nDaqFiles);
	end
end

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
	[tmpdata, t, tabs, events, info] = daqread(daqfile);
	[npts, nchan] = size(tmpdata);
	fprintf('Read %d points from file %s\n', npts, daqfile);
	if isempty(tmpdata)
		% file is empty
		fprintf('File %s is empty, terminating data loading\n\n', daqfile);
		if strcmpi(calmode, 'tones')
			calfreqs = calfreqs(1:(n-1));
		end
		break
	end	
	% ASSUME (!!) that stimulus data are collected on channel 1 (NI 0)
	% and mic data are on channel 2 (NI AI0)
	stimdata = tmpdata(:, 1);
	micdata = tmpdata(:, 2);
	clear tmpdata
	% get sample rate from the DAQ info object
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
	
% MjR: for measuring RMS of narrow bandwidth noise, filter around the CF of the
% noise prior to computing RMS below
	%------------------------------------------------------------------------
	% get highpass and lowpass filters for processing the data
	%------------------------------------------------------------------------
	% filter coefficients
	[fcoeffbHI, fcoeffaHI] = butter(forder, 3000/Fs, 'high');
	[fcoeffbLO, fcoeffaLO] = butter(forder, 5000/Fs, 'low');

	%------------------------------------------------------------------------
	% now process data
	%------------------------------------------------------------------------
	% window and filter the data
	micdataHI = filter(fcoeffbHI, fcoeffaHI, micdata);
	micdataHILO = filter(fcoeffbLO, fcoeffaLO, micdataHI);
	micdata = micdataHILO;
	
	
	%--------------------------------
	% process data according to mode
	%--------------------------------
	switch lower(calmode)
		%--------------------------------
		% WINDOW
		%--------------------------------
		case 'window'
			[rms_vals, rmw_windows] =  processWindows(micdata, rms_windowsize_ms, Fs);
			% convert to dB SPL
			dbvals{n} = dbspl(VtoPa * rms_vals);

			% decimate data for plotting
			micdata_reduced = decimate(micdata, DeciFactor);
			Fs_reduced = Fs / 10;
			% build time vectors for plotting
			t1 = ((1:length(micdata_reduced)) - 1) / Fs_reduced;
			t2 = rms_windowsize_ms * 0.001 * (0:size(dbvals{n},1)-1);
% 			t2 = rms_windowsize_ms * 0.001 * (0:rmsIndex);
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

		%--------------------------------
		% TONES
		%--------------------------------
		case 'tones'		
			[mags(n), phis(n), calfreqs(n)] = processTones(micdata, Fs, calfreqs(n), FreqDetectWidth);
			
		%--------------------------------
		% RMS
		%--------------------------------
		case 'rms'
			startidx = winstart/1000*Fs +1;
			rms_vals(n) = rms(micdata(startidx:end));
			dbvals(n) = dbspl(VtoPa * rms_vals(n));
		
		figure; 
		subplot(2,1,1)
		spectrogram(micdata(startidx:end),2^10, 2^9, 2^10, Fs, 'yaxis'); 
		ylim([0 10e4]); title('Last Measured Signal (windowed)', 'fontweight','bold')
		subplot(2,1,2)
		plot(micdata(startidx:end))
	end
	
end


% assign output vars

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
		labels = cell(ndatums, 1);
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

	case 'rms'
		out.dbvals = dbvals;
		out.rms_vals = rms_vals;
		
		figure
		subplot(211)
		plot(dbvals, '.-')
		set(gca, 'XTick', 1:length(dbvals));
		set(gca, 'XTickLabel', '');
		xlim([0 length(dbvals)+1])
		grid
		title(basename);
		ylabel('dB SPL')

		subplot(212)
		plot(rms_vals, '.-')
		set(gca, 'XTick', 1:length(dbvals));
% 		labels = cell(ndatums, 1);
% 		for n = 1:ndatums
% 			labels{n} = sprintf('%.1f', 0.001*calfreqs(n));
% 		end
% 		set(gca, 'XTickLabel', labels);
		xlim([0 length(dbvals)+1])
		grid
		ylabel('rms value');
		xlabel('The signals you presented, sequentially');
		
		% Plot spectrogram and amplitude plot of the last recorded signal
		figure
		% Assumes that signal from microphone was collected via input AI 1
		subplot(4,1,1)
		spectrogram(micdata,2^10, 2^9, 2^10, Fs, 'yaxis'); 
		ylim([0 10e4]); title('Last Recorded Signal', 'fontweight','bold')
		subplot(4,1,2)
		plot(micdata)
		% Just analyzed window
		subplot(4,1,3)
		spectrogram(micdata(startidx:end),2^10, 2^9, 2^10, Fs, 'yaxis'); 
		ylim([0 10e4]); title('Last Measured Signal (windowed)', 'fontweight','bold')
		subplot(4,1,4)
		plot(micdata(startidx:end))
% 		keyboard
		
end

out.Fs = Fs;
out.files = daqFiles;
out.path = basepath;
out.fcutoff = fcutoff;
out.forder = forder;
out.DeciFactor = DeciFactor;
out.fcoeff.a = fcoeffa;
out.fcoeff.b = fcoeffb;
out.VtoPa = VtoPa;
varargout{1} = out;

end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% processTones function
%------------------------------------------------------------------------
%------------------------------------------------------------------------
function [mag, phi, freq] = processTones(micdata, Fs, calfreq, FreqDetectWidth)
	% get spectrum of data
	[tmpfreqs, tmpmags, fmax, magmax] = daqdbfft(micdata, Fs, length(micdata));
	
	if calfreq == 0
		% if  calfreq == 0, use fmax to detect magnitude (automatic freq. detection)
		freq = fmax;
	elseif FreqDetectWidth == 0
		freq = calfreq;
	else
		% otherwise, search in a range around the provided frequency
		freqindx = find(between(tmpfreqs, calfreq - FreqDetectWidth, calfreq + FreqDetectWidth));
		[~, maxindx] = max(tmpmags(freqindx));
		freq = tmpfreqs(freqindx(maxindx));
	end
	% compute pll magnitude and phase at this frequency
	[mag, phi] = fitsinvec(micdata, 1, Fs, freq);

end
%------------------------------------------------------------------------
%------------------------------------------------------------------------

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
