%--------------------------------------------------------------------------
% NICal_NIinit_triggeredacq.m
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
% Settings/Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_Constants;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Initialize the NI device
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
iodev.NI = nidaq_ai_init('NI', iodev.Dnum);

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
% set TriggerType to HwDigital which sets triggering to
% digital type for NI hardware
set(iodev.NI.ai, 'TriggerType', 'HwDigital');
% set trigger source to the PFI0 trigger/counter input
set(iodev.NI.ai, 'HwDigitalTriggerSource', 'PFI0');
% trigger on positive-going part of signal
set(iodev.NI.ai, 'TriggerCondition', 'PositiveEdge');
% use 4 Volt trigger level
set(iodev.NI.ai, 'TriggerConditionValue', handles.TriggerLevel);
% only 1 "sweep" per trigger event 
set(iodev.NI.ai, 'TriggerRepeat', 0);
% set SamplesPerTrigger to # of samples to collect for each trigger event
set(iodev.NI.ai, 'SamplesPerTrigger', ms2samples(handles.cal.SweepDuration, iodev.Fs));

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% file logging
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set up automatic file logging
set(iodev.NI.ai, 'LogFileName', fullfile(OutputDataPath, OutputDataFile));
%-------------------------------------------------------
% set logging mode
%	'Disk'	sets logging mode to a file on disk (specified by 'LogFileName)
%	'Memory'	sets logging mode to memory only
%	'Disk&Memory'	logs to file and memory
%-------------------------------------------------------
set(iodev.NI.ai, 'LoggingMode', 'Disk&Memory');
% Disk logging mode
%	'Index' will append 00, 01, 02 etc on successive start commands
%	'Overwrite' will overwrite current log file
set(iodev.NI.ai, 'LogToDiskMode', 'Index');

TDTINIT = 1;


