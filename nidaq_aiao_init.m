function NI = nidaq_aiao_init(iface, Dnum)
%--------------------------------------------------------------------------
% NI = nidaq_aiao_init.m
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% initializes nidaq system for analog input and output (2 channels of each)
%--------------------------------------------------------------------------
% Input Arguments:
% 	iface		must be:
% 				'NI' (for traditional interface) 
% 					-OR-
% 				'NI-SESSSION' (for new, session interface)
% 					
%	Dnum		device id (usually 'Dev1')
% 
% Output Arguments:
% 	NI		struct containing settings for requested type
% 		NI.ao		analog output object
% 		NI.ai		analog input object
% 		NI.chO	analog output channel object
% 		NI.chI	analog input channel object
%--------------------------------------------------------------------------
% See also: NICal
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
%	19 Jul 2012 (JS): named to nidaq_aiao_init
%	18 Jan 2017 (SJS): updated comments
%	1 Feb 2017 (SJS): beginning transition to session DAQ interface
%--------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% some checks
%------------------------------------------------------------------------
%------------------------------------------------------------------------
disp('...starting NI hardware...');
if ~any(strcmpi(iface, {'NI', 'NI_SESSION'}))
	errordlg('%s: invalid interface %s', mfilename, iface);
	error('unknown interface');
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Now, Initialize the NI board (PCIe-6351)
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Legacy Interface
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if strcmpi(iface, 'NI')
	%---------------------------------------------------------------------
	% 'nidaq' specifies the national instruments device with traditional
	% DAQ Toolbox interface, Device number 1 (get this from the 
	% NI Measurement & Automation Explorer (a.k.a., MAX) program)
	%---------------------------------------------------------------------
	%---------------------------------------------------------------------
	% CONFIGURE ANALOG INPUT SUBSYSTEM
	%---------------------------------------------------------------------
	fprintf('Initializing NIDAQ device for analog input...')
	try
		ai = analoginput('nidaq', Dnum);
		fprintf('...done\n')
	catch errEvent
		fprintf('\nProblem while initializing NIDAQ device!\n\n')
		disp(errEvent)
		return
	end
	% create AI channels
	fprintf('creating analog input channel...')
	chI = addchannel(ai, [0 1]);
	fprintf('...done\n');
	ai.Channel(1).ChannelName = 'responseL';
	ai.Channel(2).ChannelName = 'responseR';
	%---------------------------------------------------------------------
	% CONFIGURE ANALOG OUTPUT SUBSYSTEM
	%---------------------------------------------------------------------
	fprintf('Initializing NIDAQ device for analog output...')
	try
		ao = analogoutput('nidaq', 'Dev1');
		fprintf('...done\n')
	catch errEvent
		fprintf('\nProblem while initializing NIDAQ device for output!\n\n')
		disp(errEvent)
		return
	end
	% create AO channel
	fprintf('creating analog output channels...')
	chO = addchannel(ao, [0 1]);
	fprintf('...done\n');
	ao.Channel(1).ChannelName = 'stimulusL';
	ao.Channel(2).ChannelName = 'stimulusR';
	%---------------------------------------------------------------------
	% set logging mode
	%	'Disk'	sets logging mode to a file on disk 
	%					(specified by 'LogFileName)
	%	'Memory'	sets logging mode to memory only
	%	'Disk&Memory'	logs to file and memory
	%---------------------------------------------------------------------
	set(ai, 'LoggingMode', 'Memory');
	%---------------------------------------------------------------------
	% save in NI struct and return
	%---------------------------------------------------------------------
	NI.ai = ai;
	NI.ao = ao;
	NI.chI = chI;
	NI.chO = chO;
	NI.status = 1;
	return
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Session Interface
%------------------------------------------------------------------------
%------------------------------------------------------------------------
elseif strcmpi(iface, 'NI-SESSION')
	%---------------------------------------------------------------------
	% CREATE SESSSION
	%---------------------------------------------------------------------
	% 'ni' specifies the national instruments device session
	% Device ID can be found using the daq.getDevices command
	%---------------------------------------------------------------------
	fprintf('Creating session interface...')
	try
		NI.S = daq.createSession('NI');
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
		NI.chI(1) = addAnalogInputChannel(NI.S, Dnum, 0, 'Voltage');
		NI.chI(2) = addAnalogInputChannel(NI.S, Dnum, 1, 'Voltage');
		fprintf('...done\n');
	catch errEvent
		fprintf('\nProblem while adding ai channels to NIDAQ device!\n\n')
		disp(errEvent)
		return
	end
	%---------------------------------------------------------------------
	% CONFIGURE ANALOG OUTPUT SUBSYSTEM
	%---------------------------------------------------------------------
	fprintf('Adding analog output channels...')
	try
		NI.chO(1) = addAnalogOutputChannel(NI.S, Dnum, 0, 'Voltage');
		NI.chO(2) = addAnalogOutputChannel(NI.S, Dnum, 1, 'Voltage');
		fprintf('...done\n');
	catch errEvent
		fprintf('\nProblem while adding ao channels to NIDAQ devicet!\n\n')
		disp(errEvent)
		return
	end
	% set outputs to 0
	outputSingleScan(NI.S, [0 0]);
	%---------------------------------------------------------------------
	% save in NI struct
	%---------------------------------------------------------------------
else
	errordlg('Unknown interface')
	error('Unknown interface');
end

