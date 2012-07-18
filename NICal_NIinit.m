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
%--------------------------------------------------------------------------

disp('...starting NI hardware...');

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Initialize the NI device
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%iodev.NI = handles.initfunction('NI', iodev.Dnum);
iodev.NI = nidaq_aiao_init('NI', iodev.Dnum);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set sample rate to value specified in cal settings
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%------------------------------------------------------
% AI subsystem
%------------------------------------------------------
set(iodev.NI.ai, 'SampleRate', handles.cal.Fs);
ActualRate = get(iodev.NI.ai, 'SampleRate');
if handles.cal.Fs ~= ActualRate
	warning('NICal:NIDAQ', 'Requested ai Fs (%f) ~= ActualRate (%f)', handles.cal.Fs, ActualRate);
end
iodev.Fs = ActualRate;
handles.cal.Fs = ActualRate;

%------------------------------------------------------
% AO subsystem
%------------------------------------------------------
set(iodev.NI.ao, 'SampleRate', iodev.Fs);
ActualRate = get(iodev.NI.ao, 'SampleRate');
if iodev.Fs ~= ActualRate
	warning('NICAl:NIDAQ', 'ao: Requested SampleRate (%f) ~= ActualRate (%f)', iodev.Fs, ActualRate);
end
iodev.Fs = ActualRate;
handles.cal.Fs = ActualRate;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set input range
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
for n = 1:length(iodev.NI.ai.Channel)
	iodev.NI.ai.Channel(n).InputRange = [-5 5];
end

%------------------------------------------------------------------------
% EVENT and CALLBACK PARAMETERS
%------------------------------------------------------------------------
% % first, set the object to call the SamplesAcquiredFunction when
% % BufferSize # of points are available
% set(iodev.NI.ai, 'SamplesAcquiredFcnCount', AcqSamples);
% % provide callback function handle (ai_plotpeek_callback.m)
% set(iodev.NI.ai, 'SamplesAcquiredFcn', {@ai_plotpeek_2chan_callback});

%------------------------------------------------------------------------
% HARDWARE TRIGGERING
%------------------------------------------------------------------------
% set TriggerType to manual (to synchronize ai and ao)
set([iodev.NI.ai iodev.NI.ao], 'TriggerType', 'Manual');
% set manual trigger to HW on
set(iodev.NI.ai,'ManualTriggerHwOn','Trigger')
% only 1 "sweep" per trigger event 
set(iodev.NI.ai, 'TriggerRepeat', 0);
% set SamplesPerTrigger to Inf for continous acquisition or 
% to # of samples to collect for each trigger event
set(iodev.NI.ai, 'SamplesPerTrigger', ms2samples(cal.SweepDuration, iodev.Fs));

%-------------------------------------------------------
% set logging mode
%	'Disk'	sets logging mode to a file on disk (specified by 'LogFileName)
%	'Memory'	sets logging mode to memory only
%	'Disk&Memory'	logs to file and memory
%-------------------------------------------------------
set(iodev.NI.ai, 'LoggingMode', 'Memory');

TDTINIT = 1;


