% Test Triggered Acq

% snippets to acquire in milliseconds
AcqDuration = 1000;
% sample rate in samples/second
Fs = 10000;

%---------------------------------------------------------------------
% CREATE SESSION
%---------------------------------------------------------------------
% 'ni' specifies the national instruments device session
% Device ID can be found using the daq.getDevices command
%---------------------------------------------------------------------
fprintf('Creating session interface...')
try
	S = daq.createSession('NI');
	fprintf('...done\n');
catch errEvent
	fprintf('\nCould not create session\n\n');
	disp(errEvent);
	error('Could not create event');
end
%---------------------------------------------------------------------
% CONFIGURE ANALOG INPUT SUBSYSTEM
%---------------------------------------------------------------------
fprintf('Adding analog input channels...')
try
	chI(1) = addAnalogInputChannel(S, 'Dev1', 0, 'Voltage');
	chI(2) = addAnalogInputChannel(S, 'Dev1', 1, 'Voltage');
	fprintf('...done\n');
catch errEvent
	fprintf('\nProblem while adding ai channels to NIDAQ device!\n\n')
	disp(errEvent)
	return
end
%------------------------------------------------------
% AI subsystem
%------------------------------------------------------
% set sample rate
S.Rate = Fs;
% check actual rate for mismatch
ActualRate = S.Rate;
if 1000 ~= ActualRate
	warning('testtriggereddaacq:NIDAQ', 'Requested ai Fs (%f) ~= ActualRate (%f)', ...
					handles.cal.Fs, ActualRate);
end
% store actual rate
Fs = ActualRate;
% set input range
%	range needs to be in [RangeMin RangeMax] format
aiaoRange = 5 * [-1 1];
for n = 1:length(chI)
	% set analog input range
	chI(n).Range = aiaoRange;
	% set input TerminalConfig to 'SingleEnded' (default is
	% 'Differential')
	chI(n).TerminalConfig = 'SingleEnded';
end
% make session continuous/not continuous
S.IsContinuous = false;

%------------------------------------------------------------------------
% add trigger, configure
%------------------------------------------------------------------------
addTriggerConnection(S, 'External', ...
													'Dev1/PFI0', ...
													'StartTrigger');
% trigger on rising edge
S.Connections(1).TriggerCondition = 'RisingEdge';
% Trigger timeout (seconds)
S.ExternalTriggerTimeout = 5;
% set TriggersPerRun to Inf to allow unlimited number of triggers
S.TriggersPerRun = Inf;

reset(figure(1));
fH = figure(1);
aH = axes(fH);
lH = plot(aH, 1:1000, zeros(1000, 2));
ylim(aH, [-1 1]);

%------------------------------------------------------------------------
% EVENT and CALLBACK PARAMETERS
%------------------------------------------------------------------------
% add listener for DataAvailable trigger
hl = addlistener(S, 'DataAvailable', ...
											@(src, event) testcallback(src, event, lH));
% set the object to trigger the DataAvailable callback when
% BufferSize # of points are available
S.NotifyWhenDataAvailableExceeds = ms2samples(AcqDuration, Fs);
S.DurationInSeconds = 0.001 * AcqDuration;
%---------------------------------------------------------------------
%% run!
%---------------------------------------------------------------------
fprintf('Ready to acquire data...\n\n\n');
% start in background
S.startBackground()
% wait for timeout
S.wait();
fprintf('done...\n\n\n');

%-----------------------------------------------------------------------
%% Close hardware
%-----------------------------------------------------------------------
% stop acquiring
fprintf('... terminating acquisition \n\n')
% stop acquisition
S.stop();
% stop continuous acq
S.IsContinuous = false;
% delete callback
delete(hl);
% release hardware
release(S);


% clear
clear S;

