function [handles, init_status] = NICal_NIinit_triggeredacq(handles)
%--------------------------------------------------------------------------
% [handles, init_status] = NICal_NIinit_triggeredacq(handles)
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% sets up NI data acquisition toolbox parameters
%
%	This is used for situations in which only acquisition is needed and is
%	triggered by external TTL pulse
%--------------------------------------------------------------------------
% Input Arguments:
% 	handles	struct of GUI handles from NICal
% 
% Output Arguments:
% 	handles			struct of GUI handles from NICal
%	init_status		0 if unsuccessful (handles.NI will be empty)
%						1 if successful	
%--------------------------------------------------------------------------
% See also: NICal, NICal_RunTriggeredCalibration
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 9 July 2012 (SJS)
% 				Created from SpeakerCal_tdtinit.m
% 
% Revisions:
%	9 July, 2012 (SJS) renamed for NICal project
%	15 Aug 2012 (SJS) updated comments, functionalized
%	18 Jan 2017 (SJS): updated comments
%	7 Feb 2017 (SJS): updating for Session DAQ interface
%--------------------------------------------------------------------------

fprintf('%s: starting NI hardware...\n', mfilename);
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Settings/Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_Constants;
init_status = 0;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Initialize the NI device
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
try
	if handles.DAQSESSION
		handles.iodev.NI = handles.initfunction('NI-SESSION', handles.iodev.Dnum);
	else
		handles.iodev.NI = handles.initfunction('NI', handles.iodev.Dnum);
	end
catch errMsg
	errordlg('error initializing NI device')
	warning('%s: error initializing NI device', mfilename)
	fprintf('%s: %s\n%s\n', mfilename, errMsg.identifier, errMsg.message);
	init_status = 0;
	return
end

% Different behavior depending on interface

if handles.DAQSESSION
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	% SESSSION interface
	%-----------------------------------------------------------------------
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
	handles.iodev.NI.S.IsContinuous = false;
	%------------------------------------------------------------------------
	% add trigger, configure
	%------------------------------------------------------------------------
% 	addTriggerConnection(handles.iodev.NI.S, 'External', ...
% 															'Dev1/PFI0', ...
% 															'StartTrigger');
	addTriggerConnection(handles.iodev.NI.S, 'External', ...
															'Dev1/RTSI0', ...
															'StartTrigger');
	% trigger on rising edge
	handles.iodev.NI.S.Connections(1).TriggerCondition = 'RisingEdge';
	%------------------------------------------------------------------------
	% EVENT and CALLBACK PARAMETERS
	%------------------------------------------------------------------------
	% add listener for DataAvailable trigger
	handles.iodev.hl = addlistener(handles.iodev.NI.S, 'DataAvailable', ...
												@triggered_sessioncallback);
	% set the object to trigger the DataAvailable callback when
	% BufferSize # of points are available
	handles.iodev.NI.S.NotifyWhenDataAvailableExceeds = ...
					ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);
	%-----------------------------------------------------------------------
	% Trigger timeout
	%-----------------------------------------------------------------------
	handles.iodev.NI.S.ExternalTriggerTimeout = ...
										handles.cal.TriggerSettings.TriggerTimeout;
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	
%-----------------------------------------------------------------------
% end of SESSION
%-----------------------------------------------------------------------
else
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	% LEGACY
	%-----------------------------------------------------------------------
	%-----------------------------------------------------------------------
	%------------------------------------------------------
	% AI subsystem
	%------------------------------------------------------
	% set AI sample rate
	set(handles.iodev.NI.ai, 'SampleRate', handles.cal.Fs);
	% check actual rate for mismatch
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
	%------------------------------------------------------------------------
	% HARDWARE TRIGGERING
	%------------------------------------------------------------------------
	% set TriggerType to HwDigital which sets triggering to
	% digital type for NI hardware
	% set(handles.iodev.NI.ai, 'TriggerType', 'HwDigital');
	set(handles.iodev.NI.ai, 'TriggerType', ...
						handles.cal.TriggerSettings.TriggerType);
	% set trigger source to the PFI0 trigger/counter input
	% set(handles.iodev.NI.ai, 'HwDigitalTriggerSource', 'PFI0');
	switch handles.cal.TriggerSettings.TriggerType
		case 'HwDigital'
			% use digital (TTL) channel
			set(handles.iodev.NI.ai, 'HwDigitalTriggerSource', ...
					handles.cal.TriggerSettings.TriggerSource);
		case {'HwAnalogChannel', 'HwAnalogPin'}
			% use analog channel
			set(handles.iodev.NI.ai, 'TriggerChannel', ...
						handles.iodev.NI.ai.Channel( ...
									handles.cal.TriggerSettings.TriggerSource) );
		otherwise
			error('%s: invalid TriggerType %s', ...
						mfilename, handles.cal.TriggerSettings.TriggerType);
	end
	% trigger on positive-going part of signal
	set(handles.iodev.NI.ai, 'TriggerCondition', ...
							handles.cal.TriggerSettings.TriggerCondition);
	% use 4 Volt trigger level
	set(handles.iodev.NI.ai, 'TriggerConditionValue', ...
							handles.cal.TriggerSettings.TriggerLevel);
	% only 1 "sweep" per trigger event 
	set(handles.iodev.NI.ai, 'TriggerRepeat', 0);
	% set SamplesPerTrigger to # of samples to collect for each trigger event
	set(handles.iodev.NI.ai, 'SamplesPerTrigger', ...
								ms2samples(handles.cal.SweepDuration, handles.iodev.Fs));
	%------------------------------------------------------------------------
	% EVENT and CALLBACK PARAMETERS
	%------------------------------------------------------------------------
	% first, set the object to call the SamplesAcquiredFunction when
	% BufferSize # of points are available
	set(handles.iodev.NI.ai, 'SamplesAcquiredFcnCount', ...
							ms2samples(handles.cal.SweepDuration, handles.iodev.Fs));
	% provide callback function handle (ai_plotpeek_callback.m)
	set(handles.iodev.NI.ai, 'SamplesAcquiredFcn', {@plot_callback});
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% file logging
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% set up automatic file logging
	set(handles.iodev.NI.ai, 'LogFileName', ...
							fullfile(handles.OutputDataPath, handles.OutputDataFile));
	%-------------------------------------------------------
	% set logging mode
	%	'Disk'	sets logging mode to a file on disk (specified by 'LogFileName)
	%	'Memory'	sets logging mode to memory only
	%	'Disk&Memory'	logs to file and memory
	%-------------------------------------------------------
	set(handles.iodev.NI.ai, 'LoggingMode', 'Disk&Memory');
	%-------------------------------------------------------
	% Disk logging mode
	%	'Index' will append 00, 01, 02 etc on successive start commands
	%	'Overwrite' will overwrite current log file
	%-------------------------------------------------------
	set(handles.iodev.NI.ai, 'LogToDiskMode', 'Index');
	%-------------------------------------------------------
	% set channel skew mode to Equisample
	%-------------------------------------------------------
	set(handles.iodev.NI.ai, 'ChannelSkewMode', 'Equisample');
	%set(handles.iodev.NI.ai, 'ChannelSkewMode', 'Minimum');
	%set(handles.iodev.NI.ai, 'ChannelSkewMode', 'Manual');
	%set(handles.iodev.NI.ai, 'ChannelSkew', 3.0e-06);
	%-------------------------------------------------------
	% set init_status to 1
	%-------------------------------------------------------
	init_status = 1;
end
