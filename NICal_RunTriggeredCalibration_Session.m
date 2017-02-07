%--------------------------------------------------------------------------
% NICal_RunTriggeredCalibration_Session.m
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% Runs the TTL triggered speaker calibration
%--------------------------------------------------------------------------
% See also: NICal
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	7 February, 2017 NICal_RunTriggeredCalibration
%
% Revisions:
%--------------------------------------------------------------------------
% To do:  
% 	-	some way to avoid global variables?
%--------------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_Constants;
%-----------------------------------------------------------------------
% global vars for access by DAQ callback function - needed to plot
%-----------------------------------------------------------------------
global	VtoPa Gain fcoeffa fcoeffb ...
			tvec_acq fvec Lacq Racq Lfft Rfft H SweepPoints side
% assign plot decimation factor to deciFactor
deciFactor = handles.cal.deciFactor;
% channel being calibrated: 1 = L, 2 = R, 3 = Both
side = handles.cal.side;
%---------------------------------------------
% Local
%---------------------------------------------
% set the COMPLETE flag to 0
COMPLETE = 0;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Initialization & init Scripts
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%---------------------------------------------
%---------------------------------------------
% Load Microphone calibration data
%---------------------------------------------
%---------------------------------------------
if read_ui_val(handles.FRenableCtrl) == 0
	DAscale = read_ui_str(handles.DAscaleCtrl, 'n');
	handles.cal.mic_fr = [];
	handles.cal.DAscale = DAscale;
else
	load(handles.cal.mic_fr_file, 'frdata');
	if ~isfield(frdata, 'DAscale')
		frdata.DAscale = frdata.calsettings.DAscale;
	end
	handles.cal.mic_fr = frdata;
end
% save the GUI handle information
guidata(hObject, handles);
%---------------------------------------------
%---------------------------------------------
% Microphone settings
%---------------------------------------------
%---------------------------------------------
% fix # gain values if Nchannels doesn't match # of gain values
if handles.cal.Nchannels  ~= length(handles.cal.MicGain)
	handles.cal.MicGain = handles.cal.MicGain(1) .* ...
										ones(1, handles.cal.Nchannels);
	update_ui_str(handles.MicGainCtrl, handles.cal.MicGain);
	guidata(hObject, handles);
end
% read in the gain on the mic preamp
Gain_dB = handles.cal.MicGain;
% convert dB to linear scale
Gain = invdb(Gain_dB);
% this is the sensitivity of the calibration mic in V / Pa
% if FR file is used, get sens. from there, 
if read_ui_val(handles.FRenableCtrl) == 1
	CalMic_sense = frdata.calsettings.CalMic_sense;
else
	% otherwise, use cal information
	CalMic_sense = handles.cal.MicSensitivity;
end
% pre-compute the V -> Pa conversion factor
VtoPa = (CalMic_sense^-1);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% setup output files - automagic triggering will write data to
% binary output file
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
[tmppath, tmpname] = fileparts(handles.cal.calfile);
tmpname = sprintf('%s.bin', tmpname);
tmpfilename = fullfile(tmppath, tmpname);
[filename, pathname] = uiputfile(	{'*.bin', 'Binary file (*.bin)'}, ...
												'Save Triggered data', ...
												tmpfilename );
if isequal(filename, 0) || isequal(pathname, 0)
	disp('Cancelling')
	return
else
	handles.OutputDataPath = pathname;
	handles.OutputDataFile = filename;
	[~, filebase] = fileparts(filename);
	handles.OutputMatFile = [filebase '.mat'];
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Start DAQ things
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
[handles, init_status] = NICal_NIinit_triggeredacq(handles);
guidata(hObject, handles);
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Define a bandpass filter for processing the data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Nyquist frequency
fnyq = handles.iodev.Fs/2;
% passband definition
handles.cal.fband = [handles.cal.InputHPFc handles.cal.InputLPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[handles.cal.fcoeffb, handles.cal.fcoeffa] = ...
					butter(handles.cal.forder, handles.cal.fband, 'bandpass');
fcoeffb = handles.cal.fcoeffb;
fcoeffa = handles.cal.fcoeffa;
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% write file header (FIGURE THIS OUT!!!!!!!!)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
fp = fopen(fullfile(handles.OutputDataPath, handles.OutputDataFile), ...
					'w');
fclose(fp);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Setup Plots
%-----------------------------------------------------------------------
% to speed up plotting, the vectors Lacq, Racq, tvec_acq, L/Rfft, fvec
% are pre-allocated and then those arrys are used as XDataSource and
% YDataSource for the respective plots
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-------------------------------------------------------
% create arrays for plotting (with null value)
%	create null response and time vector for plots, set up plots
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% fake acquired data
zeroacq = syn_null(handles.cal.SweepDuration, handles.iodev.Fs, 0);
zeroacq = downsample(zeroacq, handles.cal.deciFactor);
acqpts = length(zeroacq);
% time vector for stimulus plots
tvec_acq = 1000*(1/handles.iodev.Fs)*(0:(acqpts-1));
% compute # of points per sweep
SweepPoints = ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);
% acq
Lacq = zeroacq;
Racq = zeroacq;
% FFT
nfft = length(start_bin:end_bin);
tmp = zeros(1, nfft);
[fvec, Lfft] = daqdbfft(tmp, handles.iodev.Fs, nfft);
Rfft = Lfft;
% convert fvec to kHz
fvec = 0.001 * fvec;
clear tmp
%-------------------------------------------------------
% plot null data, save handles for time-domain plots
%-------------------------------------------------------
% stimulus
H.Lstim = plot(0);
H.Rstim = plot(0);
% response
H.Lacq = plot(handles.Lmicplot, tvec_acq, Lacq, 'g');
set(H.Lacq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Lacq');
xlabel(handles.Lmicplot, 'Time (ms)')
ylabel(handles.Lmicplot, 'V')
H.Racq = plot(handles.Rmicplot, tvec_acq, Racq, 'r');
set(H.Racq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Racq');
xlabel(handles.Rmicplot, 'Time (ms)')
%-------------------------------------------------------
% plot null data, save handles for frequency-domain plots
%-------------------------------------------------------
H.Lfft = plot(handles.Lfftplot, fvec, Lfft);
set(H.Lfft, 'XDataSource', 'fvec', 'YDataSource', 'Lfft');
xlabel(handles.Lfftplot, 'Frequency (kHz)')
ylabel(handles.Lfftplot, 'dBV')
H.Rfft = plot(handles.Rfftplot, fvec, Rfft);
set(H.Rfft, 'XDataSource', 'fvec', 'YDataSource', 'Rfft');
xlabel(handles.Rfftplot, 'Frequency (kHz)')

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Acquire data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
fprintf('Ready to acquire data...\n\n\n');

%----------------------------------------------------------------
% open a menu-panel for user to start (or cancel) monitoring
% of PF1 input for triggering of acquisition
%----------------------------------------------------------------
userResp = menu('Acquiring Triggered Data', 'Start', 'Cancel');

%----------------------------------------------------------------
% if user hit start, begin acquisition at triggers
%----------------------------------------------------------------
if userResp ~= 2
	runFLAG = 1;
	while runFLAG
		%START ACQUIRING
		start(handles.iodev.NI.ai);
		guidata(hObject, handles);
		fprintf('Writing to file %s \n', handles.iodev.NI.ai.LogFileName);
		% wait for TriggerTimeout seconds
		try
			wait(handles.iodev.NI.ai, handles.cal.TriggerSettings.TriggerTimeout);
		catch errEvent
			% if timeout occurs (more than handles.cal.TriggerTimeout seconds
			% elapse since prior trigger), set runFLAG to 0 to end looping
			fprintf('\n\n*** TIMEOUT!!! ***\n');
			runFLAG = 0;
		end
		% check for abort button press
		if read_ui_val(handles.AbortCtrl) == 1
			% if so, stop
			fprintf('abortion detected\n\n');
			runFLAG = 0;
		end
	end
	% stop acquiring
	fprintf('... terminating acquisition \n\n')
	stop(handles.iodev.NI.ai);
	guidata(hObject, handles);
%----------------------------------------------------------------
% otherwise, cancel
%----------------------------------------------------------------
else
	fprintf('Cancelling acquisition!\n\n');
end
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% get event log
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
EventLogAI = showdaqevents(handles.iodev.NI.ai);
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Exit gracefully (close TDT objects, etc)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Close hardware
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
disp('...closing NI devices...');
% delete and clear ai and ch0 object
delete(handles.iodev.NI.ai);
clear handles.iodev.NI.ai
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% save settings information to mat file
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
if userResp ~= 2
	cal = handles.cal;
	save(fullfile(handles.OutputDataPath, handles.OutputMatFile), ...
			'cal'			, ...
			'VtoPa'					, ...
			'EventLogAI'			, ...
			'-MAT' );
end
% save settings information to mat file
if DEBUG
	save(fullfile(pwd, 'NICal_EventLogs.mat'), ...
			'EventLogAI'			, ...
			'-MAT' );
end
COMPLETE = 1;
disp('Finished.')

