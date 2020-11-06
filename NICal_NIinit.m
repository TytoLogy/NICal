function [handles, init_status] = NICal_NIinit(handles)
%--------------------------------------------------------------------------
% NICal_NIinit.m
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% NI data acquisition toolbox parameters
% sets up nidaq system for NICal program
%------------------------------------------------------------------------
% Input Arguments:
% 	handles	struct of GUI handles from NICal
% 
% Output Arguments:
% 	handles		updated handles with initialized NI struct:
% 		for Legacy DAQ interface:
%				 		NI.ao		analog output object
%						NI.ai		analog input object
%						NI.chO	analog output channel object
%						NI.chI	analog input channel object
% 		for Session DAQ interface:
%						NI.S		session interface
%						NI.chO	analog output channel objects
%						NI.chI	analog input channel objects
% init_status	0 if unsuccessful (handles.NI will be empty)
% 					1 if successful
%------------------------------------------------------------------------
% See also: NICal, nidaq_aiao_init
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
%	1 Aug, 2012 (SJS)	added channel skew mode setting, use Equisample
% 			instead of 'Minimum" default value
%	15 Aug 2012 (SJS) updated comments, functionalized
%	18 Jan 2017 (SJS): updated comments
%	1 Feb 2017 (SJS): beginning transition to session DAQ interface
%--------------------------------------------------------------------------
	fprintf('%s: starting NI hardware...\n', mfilename);
	if handles.DAQSESSION == 0
		fprintf('%s: using legacy interface\n', mfilename);
		[handles, init_status] = legacy_init(handles);
	else
		fprintf('%s: using session interface\n', mfilename);
		[handles, init_status] = session_init(handles);
	end
	fprintf('done\n');
end

function [handles, init_status] = legacy_init(handles)
	%-----------------------------------------------------------------------
	% Settings/Constants
	%-----------------------------------------------------------------------
	NICal_Constants;
	%-----------------------------------------------------------------------
	% Initialize the NI device
	%-----------------------------------------------------------------------
	try
		handles.iodev.NI = handles.initfunction('NI', handles.iodev.Dnum);
	catch errMsg
		errordlg('error initializing NI device')
		fprintf('%s: %s\n%s\n', mfilename, errMsg.identifier, errMsg.message);
		init_status = 0;
		return
	end
	%-----------------------------------------------------------------------
	% set AI, AO sample rates to values specified in cal settings
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
	%------------------------------------------------------
	% AO subsystem
	%------------------------------------------------------
	% set AO sample rate
	set(handles.iodev.NI.ao, 'SampleRate', handles.iodev.Fs);
	% check actual rate for mismatch
	ActualRate = get(handles.iodev.NI.ao, 'SampleRate');
	if handles.iodev.Fs ~= ActualRate
		warning('NICAl:NIDAQ', 'ao: Requested SampleRate (%f) ~= ActualRate (%f)', ...
								handles.iodev.Fs, ActualRate);
	end
	% store actual rate
	handles.iodev.Fs = ActualRate;
	handles.cal.Fs = ActualRate;
	%-----------------------------------------------------------------------
	% set input range
	%-----------------------------------------------------------------------
	% range needs to be in [RangeMin RangeMax] format
	aiaoRange = 5 * [-1 1];
	% set analog input range (might be overkill to set 
	% InputRange, SensorRange and UnitsRange, but is seems to work)
	for n = 1:length(handles.iodev.NI.ai.Channel)
		handles.iodev.NI.ai.Channel(n).InputRange = aiaoRange;
		handles.iodev.NI.ai.Channel(n).SensorRange = aiaoRange;
		handles.iodev.NI.ai.Channel(n).UnitsRange = aiaoRange;
	end
	% set analog output range
	for n = 1:length(handles.iodev.NI.ao.Channel)
		handles.iodev.NI.ao.Channel(n).OutputRange = aiaoRange;
		handles.iodev.NI.ao.Channel(n).UnitsRange = aiaoRange;
	end
	%------------------------------------------------------------------------
	% HARDWARE TRIGGERING
	%------------------------------------------------------------------------
	% set TriggerType to manual (to synchronize ai and ao)
	set([handles.iodev.NI.ai handles.iodev.NI.ao], 'TriggerType', 'Manual');
	% set manual trigger to HW on
	set(handles.iodev.NI.ai, 'ManualTriggerHwOn','Trigger')
	% only 1 "sweep" per trigger event 
	set(handles.iodev.NI.ai, 'TriggerRepeat', 0);
	% set SamplesPerTrigger to Inf for continous acquisition or 
	% to # of samples to collect for each trigger event
	set(handles.iodev.NI.ai, 'SamplesPerTrigger', ...
						ms2samples(handles.cal.SweepDuration, handles.iodev.Fs));
	%-------------------------------------------------------
	% set logging mode
	%	'Disk'	sets logging mode to a file on disk (specified by 'LogFileName)
	%	'Memory'	sets logging mode to memory only
	%	'Disk&Memory'	logs to file and memory
	%-------------------------------------------------------
	set(handles.iodev.NI.ai, 'LoggingMode', 'Memory');
	%-------------------------------------------------------
	% set channel skew mode to Equisample
	%-------------------------------------------------------
	set(handles.iodev.NI.ai, 'ChannelSkewMode', 'Equisample');
	%-------------------------------------------------------
	% set init_status to 1
	%-------------------------------------------------------
	init_status = 1;
	%-------------------------------------------------------
	% Done!!!!
	%-------------------------------------------------------
end

function [handles, init_status] = session_init(handles)
	%-----------------------------------------------------------------------
	% Settings/Constants
	%-----------------------------------------------------------------------
	NICal_Constants;
	%-----------------------------------------------------------------------
	% Initialize the NI device
	%-----------------------------------------------------------------------
	try
		handles.iodev.NI = handles.initfunction('NI-SESSION', handles.iodev.Dnum);
	catch errMsg
		errordlg('error initializing NI device')
		fprintf('%s: %s\n%s\n', mfilename, errMsg.identifier, errMsg.message);
		init_status = 0;
		return
	end
	%-----------------------------------------------------------------------
	% set sample rate to value specified in cal settings
	%-----------------------------------------------------------------------
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
	%-----------------------------------------------------------------------
	% input, output channel properties
	%-----------------------------------------------------------------------
	% range needs to be in [RangeMin RangeMax] format
	aiaoRange = 5 * [-1 1];
	for n = 1:length(handles.iodev.NI.chI)
		% set analog input range
		handles.iodev.NI.chI(n).Range = aiaoRange;
		% set input TerminalConfig to 'SingleEnded' (default is
		% 'Differential')
		handles.iodev.NI.chI(n).TerminalConfig = 'SingleEnded';
	end
	% set analog output range
	for n = 1:length(handles.iodev.NI.chO)
		handles.iodev.NI.ao.chO(n).Range = aiaoRange;
	end
	%------------------------------------------------------------------------
	% HARDWARE TRIGGERING
	%------------------------------------------------------------------------
	% only 1 "sweep" per trigger event 
	handles.iodev.NI.S.TriggersPerRun = 1;
	%-------------------------------------------------------
	% set init_status to 1
	%-------------------------------------------------------
	init_status = 1;
	%-------------------------------------------------------
	% Done!!!!
	%-------------------------------------------------------
end
