% window of interest around click
% [-window_size/2...stim_delay...window_size/2]
WINDOW_SIZE = 1;

pathname = 'F:\Data2\Calibrate';
[filename, pathname] = uigetfile('*.dat','Read click data from file', pathname);

%%
if any([isempty(filename) ~exist(fullfile(pathname, filename), 'file')])
	% get filename and path
	[filename, pathname] = uigetfile(	...
													{	...
														'*.dat', 'data files (*.dat)'; ...
														'*.*', 'All Files (*.*)'	...
													}, ...
													'Pick a file' ...
												);
	% return if user pressed "Cancel" button
	if any([(filename == 0), (pathname == 0)])
		return
	end
end

datfile = fullfile(pathname, filename);

load(datfile, '-mat');

cal = C.caldata;
Fs = cal.settings.Fs;

%% plot data
SpectrumWindow = 512;

Ldata = C.raw(1, :);

figure(1);

% loop through reps
for r = 1:cal.reps
	% plot trace
	subplot(211)
	t = (1000/Fs)*((1:length(Ldata{r})) - 1);
	plot(t, Ldata{r})
	xlabel('Time (ms)')
	ylabel('V')
	title(sprintf('%s: rep %d', filename, r), 'Interpreter', 'none');
	grid

	% plot spectrogram
	subplot(212)
	[S, F, T, P] = spectrogram(	Ldata{r}, ...
											SpectrumWindow, ...
											[], ...
											SpectrumWindow, ...
											Fs	);
	P = 20*log10(P);
	P(P == -Inf) = min(min(P(P ~= -Inf)));	
	surf(1000*T, 0.001*F, P, 'edgecolor', 'none');
	xlim([min(t) max(t)])
	ylim(0.001*[0 Fs/2]);
	% set(handles.AdjSpectrumAxes, 'XTick', time_ticks)	
	view(0, 90);
	xlabel('Time (ms)')
	ylabel('Frequency (kHz)')
	colormap('gray')
	colorbar
	fprintf('press key to continue...')
	pause
	fprintf('\n');
end

%% calculations

roi_window = cal.settings.StimDelay + WINDOW_SIZE*[-0.5 0.5];
bins  = ms2bin(roi_window, Fs);
pts = bins(1):bins(2);
C.roi_window = roi_window;
C.dBval = zeros(cal.reps, 1);

figure(2)
% loop through reps
for r = 1:cal.reps
	% plot trace
	t = (1000/Fs)*pts;
	plot(t, Ldata{r}(pts))
	xlabel('Time (ms)')
	ylabel('V')
	title(sprintf('%s: rep %d', filename, r), 'Interpreter', 'none');
	grid
	
	rmsval = rms(Ldata{r}(pts));
	dBval = dbspl(cal.VtoPa * rmsval);
	
	text(cal.settings.StimDelay, max(Ldata{r}(pts)), sprintf(' dB = %.4f', dBval))
	hold on
	plot(cal.settings.StimDelay, max(Ldata{r}(pts)), '.r')
	hold off
	
	C.dBval(r) = dBval;
	
	fprintf('press key to continue...')
	pause
	fprintf('\n');
end

%% print results
fprintf('\n\n');
fprintf('Click Results\n');
fprintf('dB SPL [%.1f - %.1f] ms:\n', roi_window)
for r = 1:cal.reps
	fprintf('\t%.4f\n', C.dBval(r));
end
dBmean = mean(C.dBval);
dBstd = std(C.dBval);
fprintf('mean: %.4f dB SPL\n', dBmean);
fprintf('std: %.4f dB SPL\n', dBstd);
