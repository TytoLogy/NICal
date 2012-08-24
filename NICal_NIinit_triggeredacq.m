function [handles, init_status] = NICal_NIinit_triggeredacq(handles)
%--------------------------------------------------------------------------
% [handles, init_status] = NICal_NIinit_triggeredacq(handles)
%--------------------------------------------------------------------------
% sets up NI data acquisition toolbox parameters
%
%	This is used for situations in which only acquisition is needed and is
%	triggered by external TTL pulse
%--------------------------------------------------------------------------

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
	handles.iodev.NI = nidaq_ai_init('NI', handles.iodev.Dnum);
catch errMsg
	disp('error initializing NI device')
	init_status = 0;
	return
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set sample rate to value specified in cal settings
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%------------------------------------------------------
% AI subsystem
%------------------------------------------------------
set(handles.iodev.NI.ai, 'SampleRate', handles.cal.Fs);
ActualRate = get(handles.iodev.NI.ai, 'SampleRate');
if handles.cal.Fs ~= ActualRate
	warning('NICal:NIDAQ', 'Requested ai Fs (%f) ~= ActualRate (%f)', handles.cal.Fs, ActualRate);
end
handles.iodev.Fs = ActualRate;
handles.cal.Fs = ActualRate;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set input range
%-----------------------------------------------------------------------
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

%------------------------------------------------------------------------
% HARDWARE TRIGGERING
%------------------------------------------------------------------------
% set TriggerType to HwDigital which sets triggering to
% digital type for NI hardware
set(handles.iodev.NI.ai, 'TriggerType', 'HwDigital');
% set trigger source to the PFI0 trigger/counter input
set(handles.iodev.NI.ai, 'HwDigitalTriggerSource', 'PFI0');
% trigger on positive-going part of signal
set(handles.iodev.NI.ai, 'TriggerCondition', 'PositiveEdge');
% use 4 Volt trigger level
set(handles.iodev.NI.ai, 'TriggerConditionValue', handles.cal.TriggerLevel);
% only 1 "sweep" per trigger event 
set(handles.iodev.NI.ai, 'TriggerRepeat', 0);
% set SamplesPerTrigger to # of samples to collect for each trigger event
set(handles.iodev.NI.ai, 'SamplesPerTrigger', ms2samples(handles.cal.SweepDuration, handles.iodev.Fs));

%------------------------------------------------------------------------
% EVENT and CALLBACK PARAMETERS
%------------------------------------------------------------------------
% first, set the object to call the SamplesAcquiredFunction when
% BufferSize # of points are available
set(handles.iodev.NI.ai, 'SamplesAcquiredFcnCount', ms2samples(handles.cal.SweepDuration, handles.iodev.Fs));
% provide callback function handle (ai_plotpeek_callback.m)
set(handles.iodev.NI.ai, 'SamplesAcquiredFcn', {@plot_callback});

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% file logging
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set up automatic file logging
set(handles.iodev.NI.ai, 'LogFileName', fullfile(handles.OutputDataPath, handles.OutputDataFile));
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


