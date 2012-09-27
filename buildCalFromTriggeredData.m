function varargout = buildCalFromTriggeredData(varargin)
%------------------------------------------------------------------------
% cal = buildCalFromTriggeredData(basename, basepath, calfreqs)
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
% See also: NICal, processTriggeredData
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 20 September, 2012 (SJS) from processTriggeredData
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

AUTOFREQ = 0;
basepath = '';
basename = '';
calfreqs = [];


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
	calfreqs = varargin{3};
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
if isempty(calfreqs)	
	qVal = query_user('Auto-detect tone frequencies', 'y');
	if qVal == 0
		AUTOFREQ = 0;
		% check if user wants freqs loaded from file
		q2Val = query_user('Load tone frequencies from file', 'y')
		if q2Val == 1
			% load frequencies from file
			% open panel to get .daq file name
			[tmpfile, basepath] = uigetfile(...
					 {'*.daq', 'txt file (*.txt)'}, ...
					  'Pick a .txt file containing tone frequency list', DefaultPath);
			% check if user hit cancel (tmpfile, or basepath == 0)
			if isequal(tmpfile, 0) || isequal(basepath, 0)
				disp('Frequency file load cancelled...')
				q2Val = 0;
			end

			tmpfreqs = load(fullfile(basepath, tmpfile));
			keyboard

		end
		if q2Val == 0
			% ask user to provide frequencies
			fstr = '';
			fstr = query_uservalue('Enter frequencies, separated by spaces', '');
			calfreqs = str2num(fstr);
			fprintf('Calibration tone frequencies are:\n')
			fprintf('\t%.1f\n', calfreqs);
			fprintf('\n\n');
			clear fstr;
		end
	else
		AUTOFREQ = 1;
	end
end

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
		nDaqFiles = nDaqFiles - 1;
		break
	end	
	% ASSUME (!!) that stimulus data are collected on channel 1 (NI 0)
	% and mic data are on channel 2 (NI AI0)
	ch0data = tmpdata(:, 1);
	ch1data = tmpdata(:, 2);
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
	ch0data = sin2array(ch0data', 1, Fs);
	ch0data = filter(fcoeffb, fcoeffa, ch0data);
	ch1data = sin2array(ch1data', 1, Fs);
	ch1data = filter(fcoeffb, fcoeffa, ch1data);
	
	%--------------------------------
	% process data 
	%--------------------------------
	% get spectrum of data
	[ch0.freqs, ch0.mags, ch0fmax(n), ch0magmax(n)] = daqdbfft(ch0data, Fs, length(ch0data));
	[ch1.freqs, ch1.mags, ch1fmax(n), ch1magmax(n)] = daqdbfft(ch1data, Fs, length(ch1data));
	if AUTOFREQ
		% find frequency at peak magnitude
		calfreqs(n) = ch1fmax(n);
	end
	% compute pll magnitude and phase at this frequency
	[ch0mags(n), ch0phis(n)] = fitsinvec(ch0data, 1, Fs, calfreqs(n));
	[ch1mags(n), ch1phis(n)] = fitsinvec(ch1data, 1, Fs, calfreqs(n));
	% compute pll magnitude and phase at 2*f
	[ch0d1mags(n), ch0d1phis(n)] = fitsinvec(ch0data, 1, Fs, 2*calfreqs(n));
	[ch1d1mags(n), ch1d1phis(n)] = fitsinvec(ch1data, 1, Fs, 2*calfreqs(n));
	% compute pll magnitude and phase at 3*f
	[ch0d2mags(n), ch0d2phis(n)] = fitsinvec(ch0data, 1, Fs, 3*calfreqs(n));
	[ch1d2mags(n), ch1d2phis(n)] = fitsinvec(ch1data, 1, Fs, 3*calfreqs(n));
	
	% plot signals
	subplot(221)
	tvec = 1000 * (0:(length(ch0data)-1)) / Fs;
	plot(tvec, ch0data)
	subplot(222)
	plot(tvec, ch1data)
	
	% plot mags
	subplot(223)
	plot(ch0.freqs, ch0.mags)
	subplot(224)
	plot(ch1.freqs, ch1.mags)
	drawnow
end

% assign output vars
ndatums = length(calfreqs);
out.freqs = calfreqs;
out.mags = dbspl(VtoPa * rmssin * [ch0mags; ch1mags]);
out.phis = [ch0phis; ch1phis];
out.dist = dbspl(VtoPa * rmssin * [ch0d1mags; ch1d1mags]);
out.dist2 = dbspl(VtoPa * rmssin * [ch0d2mags; ch1d2mags]);
out.magsraw = cell(2, 1);
out.distraw = cell(2, 1);
out.distraw2 = cell(2, 1);
out.magsraw{1} = ch0mags;
out.magsraw{2} = ch1mags;
out.distraw{1} = ch0d1mags;
out.distraw{2} = ch1d1mags;
out.distraw2{1} = ch0d2mags;
out.distraw2{2} = ch1d2mags;

figure
subplot(211)
plot(0.001*calfreqs, out.mags, '.-')
set(gca, 'XTickLabel', '');
grid
title(basename);
ylabel('dB SPL')
xlim(0.001*[min(calfreqs) max(calfreqs)]);
legend('Ch0', 'Ch1')


subplot(212)
plot(0.001*calfreqs, unwrap(out.phis), '.-')
xlim(0.001*[min(calfreqs) max(calfreqs)]);
grid
ylabel('phase (rad)');
xlabel('frequency (kHz)');


out.Fs = Fs;
out.files = daqFiles;
out.path = basepath;
out.fcutoff = fcutoff;
out.forder = forder;
out.DeciFactor = DeciFactor;
out.fcoeff.a = fcoeffa;
out.fcoeff.b = fcoeffb;
out.calibration_settings = cal;
varargout{1} = out;

end


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
		rmsvals(rmsIndex) = rms(micdata(rmswindows(w-1):rmswindows(w)));
	end
end
