calpath = 'C:\Users\Calibrate\Data\ACElab\E111\20260611';
AnalysisWindow = [20 210];
AnalysisWindow = [];

D = processTriggeredBinData(...
	'inputfile', ...
		fullfile(calpath, 'tonecheck_0dBatten_80dBtgt_nicaldata.bin'), ...
	'freqfile', ...
		fullfile(calpath, 'Freqlist_4000-500-101000.txt'), ...
	'analysiswindow', AnalysisWindow);


%% plot response



%% plot raw data

figure(4)
for n = 1:length(D.data)
	% plot(timevec(length(D.data{n}), D.Fs, 'ms'),...
	% 		filtfilt(D.fcoeff.b, D.fcoeff.a, D.data{n}));
	if isempty(AnalysisWindow)
		dv = 1:length(D.data{n});
	else
		AnalysisBins = ms2bin(AnalysisWindow, D.Fs);
		dv = AnalysisBins(1):AnalysisBins(2);
	end
	tmpD = filtfilt(D.fcoeff.b, D.fcoeff.a, D.data{n});
	fftplot(tmpD(dv), D.Fs, figure(4));
end
%%

