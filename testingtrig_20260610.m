calpath = 'C:\Users\Calibrate\Data\ACElab\E111\20260610';
AnalysisWindow = [20 210];

D = processTriggeredBinData(...
	'inputfile', ...
		fullfile(calpath, '2026_iCONout_4-101kHz_0dBatten_nicaldata.bin'), ...
	'freqfile', ...
		fullfile(calpath, 'Freqlist_4000-500-101000.txt'), ...
	'analysiswindow', AnalysisWindow);

%%

figure(4)
for n = 1:length(D.data)
	% plot(timevec(length(D.data{n}), D.Fs, 'ms'),...
	% 		filtfilt(D.fcoeff.b, D.fcoeff.a, D.data{n}));

	AnalysisBins = ms2bin(AnalysisWindow, D.Fs);
	dv = AnalysisBins(1):AnalysisBins(2);
	tmpD = filtfilt(D.fcoeff.b, D.fcoeff.a, D.data{n});
	fftplot(tmpD(dv), D.Fs, figure(4));
end
