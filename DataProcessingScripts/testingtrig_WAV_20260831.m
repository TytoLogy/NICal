% calpath = 'C:\Users\Calibrate\Data\ACElab\E111\20260611';
% calfile = 'WAVtest_0dBatten_3sISI_nicaldata.bin';

calpath = 'F:\Data2\ACElab\E111\20260831';
calfile = '20260831_iCON_WAVstim_12dBatten.bin';

AnalysisWindow = [20 210];
AnalysisWindow = [];

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

rmswin = 2;

D = readBinData(fullfile(calpath, calfile));

[nSweeps, nChannels] = size(D.data);
fprintf('%s: read %d sweeps from %s\n', mfilename, nSweeps, inputfile);
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


%% plot response

for n = 1:length(D.data)
	% decimate data for plotting
	tmp = filtfilt(D.fcoeff.b, D.fcoeff.a, D.data{n});
	micdata_reduced = decimate(tmp, DeciFactor);
	Fs_reduced = D.Fs / DeciFactor;
	% build time vectors for plotting
	t1 = 1000*((1:length(micdata_reduced)) - 1) / Fs_reduced;
	t2 = D.rms_windowsize_ms * ((1:length(D.dbvals{n})) - 1);
	% plot!
	figure(n)
	subplot(211)
	plot(t1, micdata_reduced);
	grid
	ylabel('Volts');
	title(sprintf('Stimulus %d', n));
	subplot(212)
	plot(t2, D.dbvals{n}, 'Marker', '.', 'Color', 'r');
	ylabel('dB SPL')
	xlabel('Time (msec)');
	grid
	title(sprintf('Peak dB SPL = %.2f', max(D.dbvals{n})), ...
						'Interpreter', 'none')
	hold on
		plot(D.rms_windowsize_ms*D.peakdb(n, 2), D.peakdb(n, 1), 'go');
	hold off
	ylim([0 85]);
end

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

