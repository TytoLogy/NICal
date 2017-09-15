function figH = PlotSweepData(varargin)

if isempty(varargin)
	H.SweepFile = 'D:\Calibrate\CurrentSweep.mat';
else
	H = varargin{1};
end

if ~exist(H.SweepFile, 'file')
	errordlg(['No Sweep file found: ' H.SweepFile]);
	return
else
	load(H.SweepFile);
end

figH = figure;

subplot(221)
plot(tvec_acq, Lacq);
title('Left', 'Color', 'g');
xlabel('Time (ms)');
ylabel('Volts');
grid('on');

subplot(222)
plot(tvec_acq, Racq);
title('Right', 'Color', 'r');
xlabel('Time (ms)');
grid('on');

subplot(223)
plot(fvec, Lfft);
xlabel('Frequency (kHz)');
grid('on');

subplot(224)
plot(fvec, Rfft);
xlabel('Frequency (kHz)');
grid('on');
