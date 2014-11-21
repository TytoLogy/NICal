function varargout = processContinuousData(varargin)
%------------------------------------------------------------------------
% dbvals = processContinuousData(<mat file name (include path!)>)
%------------------------------------------------------------------------
% TytoLogy:NICal program
%------------------------------------------------------------------------
% 
% If.mat filename is provided, program will load that file.  Otherwise, dialog
% will be opened to select .mat file
%------------------------------------------------------------------------
% See also: NICal 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 7 November, 2014 (SJS)
%
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Variables and Constants Declarations
%------------------------------------------------------------------------
%------------------------------------------------------------------------
NICal_Constants;

% window size (in milliseconds) for computing rms (and dB) values
%	use smaller values for greater resolution, larger for coarse resolution
rms_windowsize_ms = .1;

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
calmode = 'tones';
resp = {};
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% parse inputs
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% loop through arguments
if nargin
	[basepath, basename, baseext] = fileparts(varargin{1});

else
	%----------------------------------
	% output data path and file
	%----------------------------------
	% DefaultPath = 'C:\Users\Calibrate\Matlab\Calibration\Data';
	DefaultPath = pwd;
	
	%----------------------------------
	% open panel to get .daq file name
	%----------------------------------
	[tmpfile, tmppath] = uigetfile(...
			 {'*.mat', 'mat output files (*.mat)'}, ...
			  'Pick a .mat file', DefaultPath);
	% check if user hit cancel (tmpfile, or basepath == 0)
	if isequal(tmpfile, 0) || isequal(basepath, 0)
		disp('Cancelled...')
		return
	else
		[basepath, basename, baseext] = fileparts(fullfile(tmppath, tmpfile));
	end	
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% find data files
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%-----------------------------------------------
% load .mat file for this calibration session
%-----------------------------------------------
load(fullfile(basepath, [basename '.mat']), '-MAT');


%-----------------------------------------------
% find the dat files that match the basename
%-----------------------------------------------
datname = [basename '.dat'];
if exist(fullfile(basepath, datname), 'file')
	% if found, read in rawdata
	fp = fopen(fullfile(basepath, datname), 'r');
	rawCal = readStruct(fp);
	rawResp = readCell(fp);
	if rawCal.CollectBackground
		rawBG = readCell(fp);
	end
	fclose(fp);
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% PROCESS DATA
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% cal information	
CalMic_sense = cal.MicSensitivity;
% pre-compute the V -> Pa conversion factor
VtoPa = (CalMic_sense^-1);
% precompute the volts -> RMS conversion factor for sinusoids (0.7071)
RMSsin = 1/sqrt(2);

%------------------------------------------------------------------------
% get a highpass filter for processing the  data
%------------------------------------------------------------------------
Fs = cal.Fs;
% Nyquist frequency
fnyq = Fs/2;
% filter coefficients
[fcoeffb, fcoeffa] = butter(forder, fcutoff/fnyq, 'high');

%------------------------------------------------------------------------
% now process data
%------------------------------------------------------------------------
% window and filter the data
micdata = sin2array(resp{1}, 1, Fs);
micdata = filter(fcoeffb, fcoeffa, resp{1});

[rms_vals, rmw_windows] = processWindows(micdata, rms_windowsize_ms, Fs);
% convert to dB SPL
dbvals = dbspl(VtoPa * rms_vals);

% decimate data for plotting
micdata_reduced = decimate(micdata, DeciFactor);
Fs_reduced = Fs / 10;
% build time vectors for plotting
t1 = ((1:length(micdata_reduced)) - 1) / Fs_reduced;
t2 = rms_windowsize_ms * 0.001 * (0:size(dbvals,1)-1);
% 			t2 = rms_windowsize_ms * 0.001 * (0:rmsIndex);
% 			t2 = rms_windowsize_ms * 0.001 * (0:rmsIndex);
t2 = rms_windowsize_ms * 0.001 * rmw_windows;
% plot!

subplot(211)
plot(t1, micdata_reduced);
grid
ylabel('Volts');
subplot(212)
plot(t2, dbvals, 'Marker', '.', 'Color', 'r');
ylabel('dB SPL')
xlabel('Time (seconds)');
grid

% %--------------------------------
% % TONES
% %--------------------------------
% [mags(n), phis(n), calfreqs(n)] = processTones(micdata, Fs, calfreqs(n), FreqDetectWidth);
% 
% %--------------------------------
% % RMS
% %--------------------------------
% rms_vals(n) = rms(micdata);
% dbvals(n) = dbspl(VtoPa * rms_vals(n));
		

		
% % assign output vars
% switch lower(calmode)
% 	case 'tones'
% 		ndatums = length(calfreqs);
% 		out.freqs = calfreqs;
% 		out.mags = mags;
% 		out.phis = phis;
% 		out.dbvals = dbspl(VtoPa * rmssin * out.mags);
% 
% 		figure
% 		subplot(211)
% 		plot(dbspl(VtoPa * rmssin * out.mags), '.-')
% 		set(gca, 'XTick', 1:ndatums);
% 		set(gca, 'XTickLabel', '');
% 		xlim([0 ndatums+1])
% 		grid
% 		title(basename);
% 		ylabel('dB SPL')
% 
% 		subplot(212)
% 		plot(unwrap(out.phis), '.-')
% 		set(gca, 'XTick', 1:ndatums);
% 		labels = cell(ndatums, 1);
% 		for n = 1:ndatums
% 			labels{n} = sprintf('%.1f', 0.001*calfreqs(n));
% 		end
% 		set(gca, 'XTickLabel', labels);
% 		xlim([0 ndatums+1])
% 		grid
% 		ylabel('phase (rad)');
% 		xlabel('frequency (kHz)');
% 
% 	case 'window'
% 		out.dbvals = dbvals;
% 		out.rms_windowsize_ms = rms_windowsize_ms;
% 
% 	case 'rms'
% 		out.dbvals = dbvals;
% 		out.rms_vals = rms_vals;
% end

out.Fs = Fs;
out.path = basepath;
out.fcutoff = fcutoff;
out.forder = forder;
out.DeciFactor = DeciFactor;
out.fcoeff.a = fcoeffa;
out.fcoeff.b = fcoeffb;
out.VtoPa = VtoPa;


out.dbvals = dbvals;
out.rms_windowsize_ms = rms_windowsize_ms;

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
