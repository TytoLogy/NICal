function varargout = processNICal(basename, basepath)
%------------------------------------------------------------------------
% dbvals = processDAQdB(basename)
%------------------------------------------------------------------------
% 
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
fcutoff = 500;
% filter order
forder = 3;

% Decimation factor - plotted data will be 1 / DeciFactor shorter
% and sampling rate will be Fs / DeciFactor
DeciFactor = 10;

%--------------------------------------------------------
% find data files
%--------------------------------------------------------
if nargin == 0
	% output data path and file
	% DefaultPath = 'C:\Users\Calibrate\Matlab\Calibration\Data';
	DefaultPath = pwd;
	  
	[tmpfile, basepath] = uigetfile(...
			 {'*.daq', 'DAQ output files (*.daq)'}, ...
			  'Pick a file', DefaultPath);

	if isequal(tmpfile, 0) | isequal(basepath, 0)
		disp('Cancelled...')
		return
	end
	
	[~, extfile, fext] = fileparts(tmpfile);
	basename = extfile(1:(length(extfile) - 2))
	
else
	if ~exist('basepath', 'var') || isempty(basepath)
		basepath = pwd;
	end
end
[daqFiles, daqNumbers] = find_daq_files(basename, basepath);
nDaqFiles = length(daqFiles);

% load .mat file for this calibration session
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
	if isempty(tmpdata)
		break
	end
	
	stimdata = tmpdata(:, 1);
	micdata = tmpdata(:, MicChannel+1);
	clear tmpdata
	% sample rate
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

	
	[tmpfreqs, tmpmags, fmax, magmax] = daqdbfft(micdata, Fs, length(micdata));

	freqs(n) = fmax;
	[mags(n), phis(n)] = fitsinvec(micdata, 1, Fs, fmax);
end

ndatums = length(freqs);
out.freqs = freqs;
out.mags = mags;
out.phis = phis;
out.dbvals = dbvals;
varargout{1} = out;

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
	labels{n} = sprintf('%.1f', 0.001*freqs(n));
end
set(gca, 'XTickLabel', labels);
xlim([0 ndatums+1])
grid
ylabel('phase (rad)');
xlabel('frequency (kHz)');
