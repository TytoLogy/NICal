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
%	23 Mar 2017 (SJS): fixing ignored handles.cal.InputFilter
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
global	VtoPa Gain fcoeffa fcoeffb deciFactor start_bin end_bin ...
			filtEnable tvec_acq fvec Lacq Racq Lfft Rfft H nfft side
% assign plot decimation factor to deciFactor
deciFactor = handles.cal.deciFactor;
% channel being calibrated: 1 = L, 2 = R, 3 = Both
side = handles.cal.Side;
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
% assign filter/no filter to filtEnable
filtEnable = handles.cal.InputFilter;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% write file header (FIGURE THIS OUT!!!!!!!!)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
datafile = fullfile(handles.OutputDataPath, handles.OutputDataFile);
fp = fopen(datafile,	'w');
writeStruct(fp, handles.cal, 'cal');
fclose(fp);

%-----------------------------------------------------------------------
% set the start and end bins for the calibration based on 
% front panel settings
%-----------------------------------------------------------------------

% 3 Nov 2020: modifying to be more general (don't assume ramped stimuli)
%{
OLD CODE
start_bin = ms2bin(handles.cal.StimDelay + handles.cal.StimRamp, ...
							handles.iodev.Fs);
if ~start_bin
	start_bin = 1;
end
end_bin = start_bin + ms2bin(handles.cal.StimDuration-handles.cal.StimRamp, ...
											handles.iodev.Fs);
%}
% compute startbin using stimulus delay
start_bin = ms2bin(handles.cal.StimDelay, handles.iodev.Fs);
% if start_bin is 0 of undefined, set it to 1
if ~start_bin
	start_bin = 1;
end
end_bin = start_bin + ms2bin(handles.cal.StimDuration, handles.iodev.Fs);

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
plotacq = downsample(zeroacq, handles.cal.deciFactor);
acqpts = length(zeroacq);
% time vector for stimulus plots
tvec_acq = 1000*(handles.cal.deciFactor/handles.iodev.Fs) * ...
															(0:(length(plotacq)-1));
% compute # of points per sweep
% SweepPoints = ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);
% acq
Lacq = plotacq;
Racq = plotacq;
% FFT
% temp = zeroacq(start_bin:end_bin);
temp = randn(size(zeroacq(start_bin:end_bin)));
nfft = length(temp);
% temp = zeros(1, nfft);
[fvec, Lfft] = daqdbfft(temp, handles.iodev.Fs, nfft);
Rfft = Lfft;
% convert fvec to kHz
fvec = 0.001 * fvec;
clear temp
%-------------------------------------------------------
% plot null sweep data, save handles for time-domain plots
%-------------------------------------------------------
% stimulus (null)
H.Lstim = plot(handles.Lstimplot, 0);
H.Rstim = plot(handles.Rstimplot, 0);
% response
H.Lacq = plot(handles.Lmicplot, tvec_acq, Lacq, 'g');
set(H.Lacq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Lacq');
xlabel(handles.Lmicplot, 'Time (ms)')
ylabel(handles.Lmicplot, 'V')
H.Racq = plot(handles.Rmicplot, tvec_acq, Racq, 'r');
set(H.Racq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Racq');
xlabel(handles.Rmicplot, 'Time (ms)')
%-------------------------------------------------------
% plot null fft data, save handles for frequency-domain plots
%-------------------------------------------------------
H.Lfft = plot(handles.Lfftplot, fvec, Lfft);
set(H.Lfft, 'XDataSource', 'fvec', 'YDataSource', 'Lfft');
xlabel(handles.Lfftplot, 'Frequency (kHz)')
ylabel(handles.Lfftplot, 'dBV')
H.Rfft = plot(handles.Rfftplot, fvec, Rfft);
set(H.Rfft, 'XDataSource', 'fvec', 'YDataSource', 'Rfft');
xlabel(handles.Rfftplot, 'Frequency (kHz)')
% assign handles for text displays
H.LValText = handles.LValText;
H.LSPLText = handles.LSPLText;
H.RValText = handles.RValText;
H.RSPLText = handles.RSPLText;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% EVENT and CALLBACK PARAMETERS
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% add listener for DataAvailable trigger
handles.iodev.hl = addlistener(handles.iodev.NI.S, 'DataAvailable', ...
				@(src, event) triggered_sessioncallback(src, event, datafile));
% set the object to trigger the DataAvailable callback when
% BufferSize # of points are available
handles.iodev.NI.S.NotifyWhenDataAvailableExceeds = ...
				ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);
			
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
% 	runFLAG = 1;
	% START ACQUIRING
	startBackground(handles.iodev.NI.S);
	guidata(hObject, handles);
	fprintf('Writing to file %s \n', ...
						fullfile(handles.OutputDataPath, handles.OutputDataFile));

	% wait for TriggerTimeout
	handles.iodev.NI.S.wait();
	fprintf('done...\n\n\n');
% 	runFLAG = 0;

% 	% wait for triggers
% 	while runFLAG
% 		try
% 			% wait for TriggerTimeout seconds
% 			handles.iodev.NI.S.wait(handles.cal.TriggerSettings.TriggerTimeout);
% 		catch errEvent
% 			% if timeout occurs (more than handles.cal.TriggerTimeout seconds
% 			% elapse since prior trigger), set runFLAG to 0 to end looping
% 			fprintf('\n\n*** TIMEOUT!!! ***\n');
% 			runFLAG = 0;
% 		end
% 		% check for abort button press
% 		if read_ui_val(handles.AbortCtrl) == 1
% 			% if so, stop
% 			fprintf('abortion detected\n\n');
% 			runFLAG = 0;
% 		end
% 	end
	guidata(hObject, handles);
%----------------------------------------------------------------
% otherwise, cancel
%----------------------------------------------------------------
else
	fprintf('Cancelling acquisition!\n\n');
end
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Exit gracefully (close TDT objects, etc)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
% Close hardware
%-----------------------------------------------------------------------
% stop acquiring
fprintf('... terminating acquisition \n\n')
% stop acquisition
handles.iodev.NI.S.stop();
% delete callback
delete(handles.iodev.hl);
% release hardware
release(handles.iodev.NI.S);
% clear
clear handles.iodev.NI.S;

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
			'-MAT' );
end
COMPLETE = 1;
disp('Finished.')

