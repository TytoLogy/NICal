calpath = 'F:\Data2\ACElab\E111\20260831';
calfile = 'E111_20260831_allWAVstim_12dBatten.bin';
calfile = 'E111_20260831_allWAVstim_6dBatten.bin';
inputfile = fullfile(calpath, calfile);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Variables and Constants Declarations
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% stimulus channel
S = 1;
% mic channel
M = 2;

AnalysisWindow = [];

% window size (in milliseconds) for computing rms (and dB) values
%	use smaller values for greater resolution, larger for coarse resolution
rms_windowsize_ms = 2;
rmswin = 2;

% highpass cutoff frequency (Hz)
fcutoff = 3000;
% filter order
forder = 3;

% Decimation factor for plotting - plotted data will be 1 / DeciFactor shorter
% and sampling rate will be Fs / DeciFactor
DeciFactor = 10;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% read in data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
D = readBinData(fullfile(calpath, calfile));
[nSweeps, nChannels] = size(D.data);
fprintf('%s: read %d sweeps from %s\n', mfilename, nSweeps, inputfile);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Get Filter for Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% get sample rate from the cal struct
Fs = D.cal.Fs;
%------------------------------------------------------------------------
% get a highpass filter for processing the  data
%------------------------------------------------------------------------
% Nyquist frequency
fnyq = Fs/2;
% filter coefficients
[fcoeffb, fcoeffa] = butter(forder, fcutoff/fnyq, 'high');
% SPL conversion
VtoPa = 1 ./ (D.cal.Gain(1) * invdb(D.cal.MicGain(1)) * D.cal.MicSensitivity);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Process Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

D.rms_windowsize_ms = rms_windowsize_ms;

%--------------------------------
% loop through "daqfiles" (now sweeps)
%--------------------------------
% for n = 1:nSweeps
for n = 50
	if nChannels == 2
		tmpdata = [D.data{n, 1} D.data{n, 2}];
	else
		% if only one channel, "fake" the second channel by duplicating the
		% first... not most efficient but is expedient.
		tmpdata = [D.data{n, 1} D.data{n, 1}];
	end
	% ASSUME (!!) that stimulus data are collected on channel 1 (NI 0)
	% and mic data are on channel 2 (NI AI0)
	stimdata = tmpdata(:, S); %#ok<NASGU>
	micdata = tmpdata(:, M);
	clear tmpdata

	%------------------------------------------------------------------------
	% now process data
	%------------------------------------------------------------------------
	% window and filter the data
	micdata = sin2array(micdata', 5, Fs);
	micdata = filtfilt(fcoeffb, fcoeffa, micdata);

	[rms_vals, ~] =  processWindows(micdata, rms_windowsize_ms, Fs);
	
	% throw away last rms value as it is likely to be corrupted by
	% truncation (rms_bins don't evenly divide nsamples)
	D.rms_vals{n} = rms_vals(1:(end-1));
	
	% convert to dB SPL
	D.dbvals{n} = dbspl(VtoPa * D.rms_vals{n}); %#ok<AGROW>

	% decimate data for plotting
	micdata_reduced = decimate(micdata, DeciFactor);
	Fs_reduced = Fs / DeciFactor;
	% build time vectors for plotting
	t1 = 1000*((1:length(micdata_reduced)) - 1) / Fs_reduced;
	t2 = rms_windowsize_ms * ((1:length(D.rms_vals{n})) - 1);
	% plot!
	figure(n)
	subplot(311)
	plot(t1, micdata_reduced);
	grid
	ylabel('Volts');
	title(sprintf('WAV Record %d', n));
	
	subplot(312)
	plot(t2, D.dbvals{n}, 'Marker', '.', 'Color', 'b');
	ylabel('dB SPL')
	xlabel('Time (msec)');
	grid
	[pk, pkbin] = max(D.dbvals{n});
	D.peakdB(n, 1) = pk;
	D.peakdB(n, 2) = pkbin;
	title(sprintf('Peak dB SPL = %.2f', D.peakdB(n, 1)), ...
				'Interpreter', 'none')
	hold on
		plot(t2(pkbin), pk, 'ro');
	hold off
	subplot(313)
	spectrogram(micdata, 256, 250, 256, Fs, 'yaxis')
	ylim([1 125])
	colorbar off
	colormap('gray')
	
	
end



%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot response
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% for n = 1:length(D.data)
% 	% decimate data for plotting
% 	tmp = filtfilt(D.fcoeff.b, D.fcoeff.a, D.data{n});
% 	micdata_reduced = decimate(tmp, DeciFactor);
% 	Fs_reduced = D.Fs / DeciFactor;
% 	% build time vectors for plotting
% 	t1 = 1000*((1:length(micdata_reduced)) - 1) / Fs_reduced;
% 	t2 = D.rms_windowsize_ms * ((1:length(D.dbvals{n})) - 1);
% 	% plot!
% 	figure(n)
% 	subplot(211)
% 	plot(t1, micdata_reduced);
% 	grid
% 	ylabel('Volts');
% 	title(sprintf('Stimulus %d', n));
% 	subplot(212)
% 	plot(t2, D.dbvals{n}, 'Marker', '.', 'Color', 'r');
% 	ylabel('dB SPL')
% 	xlabel('Time (msec)');
% 	grid
% 	title(sprintf('Peak dB SPL = %.2f', max(D.dbvals{n})), ...
% 						'Interpreter', 'none')
% 	hold on
% 		plot(D.rms_windowsize_ms*D.peakdb(n, 2), D.peakdb(n, 1), 'go');
% 	hold off
% 	ylim([0 85]);
% end

%% save peak db data

[~, tmp, ~] = fileparts(calfile);
dbfile = [tmp '_peakdB.txt'];
fp = fopen(fullfile(calpath, dbfile), 'wt');
fprintf(fp, 'Stimulus\tPeakdBSPL\tPeakTime_bin\tPeakTime_ms\n');
for n = 1:length(D.peakdb)
	fprintf(fp, '%d\t%f\t%d\t%.1f\n', n, D.peakdb(n, 1), D.peakdb(n, 2), ...
												D.rms_windowsize_ms*D.peakdb(n, 2));
end
fclose(fp);
%% plot raw data

% figure(4)
% for n = 1:length(D.data)
% 	% plot(timevec(length(D.data{n}), D.Fs, 'ms'),...
% 	% 		filtfilt(D.fcoeff.b, D.fcoeff.a, D.data{n}));
% 	if isempty(AnalysisWindow)
% 		dv = 1:length(D.data{n});
% 	else
% 		AnalysisBins = ms2bin(AnalysisWindow, D.Fs);
% 		dv = AnalysisBins(1):AnalysisBins(2);
% 	end
% 	tmpD = filtfilt(D.fcoeff.b, D.fcoeff.a, D.data{n});
% 	fftplot(tmpD(dv), D.Fs, figure(4));
% end
%%




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

