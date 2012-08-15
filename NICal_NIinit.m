
function [handles, init_status] = NICal_NIinit(handles)
%--------------------------------------------------------------------------
% NICal_NIinit.m
%--------------------------------------------------------------------------
% sets up NI data acquisition toolbox parameters
%
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
%	1 Aug, 2012 (SJS)	added channel skew mode setting, use Equisample
% 			instead of 'Minimum" default value
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
%handles.iodev.NI = handles.initfunction('NI', iodev.Dnum);
try
	handles.iodev.NI = nidaq_aiao_init('NI', handles.iodev.Dnum);
catch
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

%------------------------------------------------------
% AO subsystem
%------------------------------------------------------
set(handles.iodev.NI.ao, 'SampleRate', handles.iodev.Fs);
ActualRate = get(handles.iodev.NI.ao, 'SampleRate');
if handles.iodev.Fs ~= ActualRate
	warning('NICAl:NIDAQ', 'ao: Requested SampleRate (%f) ~= ActualRate (%f)', handles.iodev.Fs, ActualRate);
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
set(handles.iodev.NI.ai,'ManualTriggerHwOn','Trigger')
% only 1 "sweep" per trigger event 
set(handles.iodev.NI.ai, 'TriggerRepeat', 0);
% set SamplesPerTrigger to Inf for continous acquisition or 
% to # of samples to collect for each trigger event
set(handles.iodev.NI.ai, 'SamplesPerTrigger', ms2samples(handles.cal.SweepDuration, handles.iodev.Fs));

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
%set(handles.iodev.NI.ai, 'ChannelSkewMode', 'Minimum');
%set(handles.iodev.NI.ai, 'ChannelSkewMode', 'Manual');
%set(handles.iodev.NI.ai, 'ChannelSkew', 3.0e-06);

%-------------------------------------------------------
% set init_status to 1
%-------------------------------------------------------
init_status = 1;


