calpath = 'C:\Users\Calibrate\Data\ACElab\E111\20260610';
AnalysisWindow = [20 210];

D = processTriggeredBinData(...
	'inputfile', ...
		fullfile(calpath, '2026_iCONout_4-101kHz_0dBatten_nicaldata.bin'), ...
	'freqfile', ...
		fullfile(calpath, 'Freqlist_4000-500-101000.txt'), ...
	'analysiswindow', AnalysisWindow);

%% save data
caldatafile = '2026_iCONout_4-101kHz_0dBatten_caldata.mat';
caldata.freqs = D.freqs;
caldata.mags = D.mags;
caldata.phis = D.phis;
caldata.dbcals = D.dbvals;
caldata.datafile = fullfile(D.path, D.files);
caldata.cal = D.cal;
caldata.freqfile = D.freqfile;
caldata.AnalysisWindow = D.AnalysisWindow;
save(fullfile(calpath, caldatafile), 'caldata');

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
