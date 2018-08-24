
pathname = '/Users/sshanbhag/Work/Data/Audio/Calibration/CrestV450Amp/14Mar2017';
filename = 'CrestV450_0.1V_sweep.dat';
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
fp = fopen(datfile, 'r');

% read cal struct
cal = readStruct(fp);

% allocate data storage
if cal.Side == 1
	Ldata = cell(cal.Nreps, 1);
	Rdata = {};
elseif cal.Side == 2
	Ldata = {};
	Rdata = cell(cal.Nreps, 1);
elseif cal.Side == 3
	Ldata = cell(cal.Nreps, 1);
	Rdata = cell(cal.Nreps, 1);
end
% read L data
if (cal.Side == 1) || (cal.Side == 3)
	% read L sweeps
	for s = 1:cal.Nreps
		tmp = readCell(fp);
		Ldata{s} = tmp{1};
	end
end
% read R data
if (cal.Side == 2) || (cal.Side == 3)
	% read L sweeps
	for s = 1:cal.Nreps
		tmp = readCell(fp);
		Rdata{s} = tmp{1}; %#ok<SAGROW>
	end
end	
fclose(fp);

%% plot data
figure(1);
subplot(211)

t = (1000/cal.Fs)*((1:length(Ldata{1})) - 1);
plot(t, Ldata{1})
xlabel('Time (ms)')
ylabel('V')

SpectrumWindow = 512;

	
subplot(212)

[S, F, T, P] = spectrogram(	Ldata{1}, ...
										SpectrumWindow, ...
										[], ...
										SpectrumWindow, ...
										cal.Fs	);
P = 20*log10(P);
P(P == -Inf) = min(min(P(P ~= -Inf)));	
surf(1000*T, 0.001*F, P, 'edgecolor', 'none');
xlim([min(t) max(t)])
ylim(0.001*[cal.Fmin cal.Fmax]);
% set(handles.AdjSpectrumAxes, 'XTick', time_ticks)	
view(0, 90);
xlabel('Time (ms)')
ylabel('Frequency (kHz)')
colormap('gray')
colorbar


