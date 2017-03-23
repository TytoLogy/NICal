%------------------------------------------------------------------------
% NICal_Monitor.m
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 28 August, 2012 (SJs)
% 	- adapted from monitordB script
%
% Revisions:
%	7 Feb 2017 (SJS): updated for DAQ Session interface
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Constants
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% NICal
NICal_Constants;
% global (needed for use by callback)
global VtoPa Gain fcoeffa fcoeffb filtEnable ...
		tvec_acq fvec Lacq Racq Lfft Rfft H SweepPoints

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Need to do different things depending on state of button
%------------------------------------------------------------------------
%------------------------------------------------------------------------
currentState = read_ui_val(handles.MonitorCtrl);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%***** start monitor
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if currentState == 1
	%---------------------------------------------
	% Update Monitor Button
	%---------------------------------------------
	update_ui_str(handles.MonitorCtrl, 'Monitor ON');
	set(handles.MonitorCtrl, 'FontAngle', 'italic');
	% save the GUI handle information
	guidata(hObject, handles);
	%---------------------------------------------
	% Load Microphone calibration data
	%---------------------------------------------
	if read_ui_val(handles.FRenableCtrl) == 0
		DAscale = read_ui_str(handles.DAscaleCtrl, 'n');
		handles.cal.mic_fr = [];
		handles.cal.DAscale = DAscale;
	else
		load(handles.cal.mic_fr_file, 'frdata', '-MAT');
		if ~isfield(frdata, 'DAscale')
			frdata.DAscale = frdata.calsettings.DAscale;
		end
		handles.cal.mic_fr = frdata;
	end
	%-------------------------------------------------------------
	% need some conversion factors
	%-------------------------------------------------------------
	Gain = handles.cal.Gain;
	% this is the sensitivity of the calibration mic in V / Pa
	% if FR file is used, get sens. from there, 
	if read_ui_val(handles.FRenableCtrl) == 1
		CalMic_sense = frdata.calsettings.MicSensitivity;
	else
		% otherwise, use cal information
		CalMic_sense = handles.cal.MicSensitivity;
	end
	% pre-compute the V -> Pa conversion factor
	VtoPa = (CalMic_sense^-1);
	
	%-------------------------------------------------------------
	% Start DAQ things
	%-------------------------------------------------------------
	% Initialize the NI device
	try
		if handles.DAQSESSION
			handles.iodev.NI = nidaq_ai_init('NI-SESSION', handles.iodev.Dnum);
		else
			handles.iodev.NI = nidaq_ai_init('NI', handles.iodev.Dnum);
		end
	catch errMsg
		warning('%s: error initializing NI device', mfilename)
		fprintf('Message: %s\n', errMsg.message');
		init_status = 0;
		update_ui_str(handles.MonitorCtrl, 'Monitor');
		update_ui_val(handles.MonitorCtrl, 0);
		set(handles.MonitorCtrl, 'FontAngle', 'normal');
		return
	end
	
	%-------------------------------------------------------------
	% disable some front panel controls
	%-------------------------------------------------------------
	disable_ui(handles.RunCalibrationCtrl);
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Hardware Settings
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	if handles.DAQSESSION
		%-----------------------------------------------------------------------
		%-----------------------------------------------------------------------
		% SESSSION interface
		%-----------------------------------------------------------------------
		%------------------------------------------------------
		% AI subsystem
		%------------------------------------------------------
		% set sample rate
		handles.iodev.NI.S.Rate = handles.cal.Fs;
		% check actual rate for mismatch
		ActualRate = handles.iodev.NI.S.Rate;
		if handles.cal.Fs ~= ActualRate
			warning('NICal:NIDAQ', 'Requested ai Fs (%f) ~= ActualRate (%f)', ...
							handles.cal.Fs, ActualRate);
		end
		% store actual rate
		handles.iodev.Fs = ActualRate;
		handles.cal.Fs = ActualRate;
		% set input range
		%	range needs to be in [RangeMin RangeMax] format
		aiaoRange = 5 * [-1 1];
		for n = 1:length(handles.iodev.NI.chI)
			% set analog input range
			handles.iodev.NI.chI(n).Range = aiaoRange;
			% set input TerminalConfig to 'SingleEnded' (default is
			% 'Differential')
			handles.iodev.NI.chI(n).TerminalConfig = 'SingleEnded';
		end
		% make session continuous
		handles.iodev.NI.S.IsContinuous = true;
		%------------------------------------------------------------------------
		% EVENT and CALLBACK PARAMETERS
		%------------------------------------------------------------------------
		% add listener for DataAvailable trigger
		handles.iodev.hl = addlistener(handles.iodev.NI.S, 'DataAvailable', ...
													@monitor_sessioncallback);
		% set the object to trigger the DataAvailable callback when
		% BufferSize # of points are available
		handles.iodev.NI.S.NotifyWhenDataAvailableExceeds = ...
						ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);

	else
		%-----------------------------------------------------------------------
		%-----------------------------------------------------------------------
		% Legacy interface
		%-----------------------------------------------------------------------
		%------------------------------------------------------
		% AI subsystem
		%------------------------------------------------------
		% set sample rate
		set(handles.iodev.NI.ai, 'SampleRate', handles.cal.Fs);
		% check to see what actual rate is
		ActualRate = get(handles.iodev.NI.ai, 'SampleRate');
		if handles.cal.Fs ~= ActualRate
			warning('NICal:NIDAQ', 'Requested ai Fs (%f) ~= ActualRate (%f)', ...
									handles.cal.Fs, ActualRate);
		end
		% store actual rate
		handles.iodev.Fs = ActualRate;
		handles.cal.Fs = ActualRate;
		% set input range
		% range needs to be in [RangeMin RangeMax] format
		aiaoRange = 5 * [-1 1];
		% set analog input range (might be overkill to set 
		% InputRange, SensorRange and UnitsRange, but is seems to work)
		for n = 1:length(handles.iodev.NI.ai.Channel)
			handles.iodev.NI.ai.Channel(n).InputRange = aiaoRange;
			handles.iodev.NI.ai.Channel(n).SensorRange = aiaoRange;
			handles.iodev.NI.ai.Channel(n).UnitsRange = aiaoRange;
		end
		% set SamplesPerTrigger to Inf for continous acquisition
		set(handles.iodev.NI.ai, 'SamplesPerTrigger', Inf);
		% set TriggerType to 'Manual' so that program starts acquisition
		set(handles.iodev.NI.ai, 'TriggerType', 'Manual');
		% set input type to single ended
		set(handles.iodev.NI.ai, 'InputType', 'SingleEnded');
		%------------------------------------------------------------------------
		% EVENT and CALLBACK PARAMETERS
		%------------------------------------------------------------------------
		% first, set the object to call the SamplesAcquiredFunction when
		% BufferSize # of points are available
		set(handles.iodev.NI.ai, 'SamplesAcquiredFcnCount', ...
							ms2samples(handles.cal.SweepDuration, handles.iodev.Fs));
		% provide callback function
		set(handles.iodev.NI.ai, 'SamplesAcquiredFcn', {@monitor_callback});
		% set logging mode
		%	'Disk'			sets logging mode to a file on disk 
		%						(specified by 'LogFileName)
		%	'Memory'			sets logging mode to memory only
		%	'Disk&Memory'	logs to file and memory
		set(handles.iodev.NI.ai, 'LoggingMode', 'Memory');
		% set channel skew mode to Equisample
		set(handles.iodev.NI.ai, 'ChannelSkewMode', 'Equisample');
	end
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Plots Initialization
	% create null acq and time vectors for plots, set up plots
	%-----------------------------------------------------------------------
	% to speed up plotting, the vectors Lacq, Racq, tvec_acq, L/Rfft, fvec
	% are pre-allocated and then those arrys are used as XDataSource and
	% YDataSource for the respective plots
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	% sample interval
	dt = 1/handles.iodev.Fs;
	zeroacq = syn_null(handles.cal.SweepDuration, handles.iodev.Fs, 0);
	SweepPoints = length(zeroacq);
	% time vector for stimulus plots
	tvec_acq = 1000*(1/handles.iodev.Fs)*(0:(SweepPoints-1));
	%-------------------------------------------------------
	% create arrays for plotting and plot them
	%-------------------------------------------------------
	% acq
	Lacq = zeroacq;
	Racq = zeroacq;
	% FFT
	[fvec, Lfft] = daqdbfft(zeroacq, handles.iodev.Fs, SweepPoints);
	Rfft = Lfft;
	% convert fvec to kHz
	fvec = 0.001 * fvec;
	clear tmp
	%----------------------------------------------------------------
	% plot null data, save handles in H struct for time-domain plots
	%----------------------------------------------------------------
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
	xlabel(handles.Rfftplot, 'Frequency (kHz)');
	%-------------------------------------------------------
	% store text handles in H
	%-------------------------------------------------------
	H.LValText = handles.LValText;
	H.LSPLText = handles.LSPLText;
	H.RValText = handles.RValText;
	H.RSPLText = handles.RSPLText;
	%-------------------------------------------------------
	% update handles
	%-------------------------------------------------------
	guidata(hObject, handles);
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define a bandpass filter for processing the data
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Nyquist frequency
	fnyq = handles.iodev.Fs / 2;
	% passband definition
	fband = [handles.cal.InputHPFc handles.cal.InputLPFc] ./ fnyq;
	% filter coefficients using a butterworth bandpass filter
	[fcoeffb, fcoeffa] = ...
						butter(handles.cal.forder, fband, 'bandpass');
	% assign filter/no filter to filtEnable
	filtEnable = handles.cal.InputFilter;
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% START ACQUIRING
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	if handles.DAQSESSION
		% start in background
		fprintf('Starting Monitor\n');
		startBackground(handles.iodev.NI.S);
	else
		start(handles.iodev.NI.ai);
		trigger(handles.iodev.NI.ai);
	end
	guidata(hObject, handles);
	
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%***** stop monitor
%------------------------------------------------------------------------
%------------------------------------------------------------------------
else
	update_ui_str(handles.MonitorCtrl, 'Monitor');
	set(handles.MonitorCtrl, 'FontAngle', 'normal');
	guidata(hObject, handles);
	%------------------------------------------------------------------------
	% clean up
	%------------------------------------------------------------------------
	disp('...closing NI devices...');
	if handles.DAQSESSION
		fprintf('Stopping Monitor\n');
		% stop acquisition
		handles.iodev.NI.S.stop();
		% stop continuous acq
		handles.iodev.NI.S.IsContinuous = false;
		% delete callback
		delete(handles.iodev.hl);
		% release hardware
		release(handles.iodev.NI.S);
		% clear
		clear handles.iodev.NI.S;
	else
		% stop acquiring
		stop(handles.iodev.NI.ai);
		% delete and clear ai and ch0 object
		delete(handles.iodev.NI.ai);
		clear handles.iodev.NI.ai
	end
	
	% re-enable the RunCalibration button
	enable_ui(handles.RunCalibrationCtrl);
end	% END if currentState
