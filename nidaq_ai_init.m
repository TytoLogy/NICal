function NI = nidaq_ai_init(iface, Dnum)
%--------------------------------------------------------------------------
% NI = nidaq_ai_init.m
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% initializes nidaq system for analog acq
% designed for use with National Instruments' PCIe-6351 board
%------------------------------------------------------------------------
% Input Arguments:
% 	iface		must be:
% 				'NI' (for traditional interface) 
% 					-OR-
% 				'NI-SESSSION' (for new, session interface)
%
%	Dnum		device id (usually 'Dev1')
% 
% Output Arguments:
% 	NI		struct containing interface for hardware
%		SESSION:
%			NI.S		DAQ Toolbox Session interface object
%			NI.chI	analog input channel object
%		LEGACY:
%			NI.ai		analog input object (LEGACY only)
%			NI.chI	analog input channel object
%------------------------------------------------------------------------
% See also: NICal, nidaq_aiao_init, DAQ Toolbox (Matlab)
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 18 July 2012 (SJS)
% 				Created from nidaq_init.m
% 
% Revisions:
%	18 Jan 2017 (SJS): updated comments
%	1 Feb 2017 (SJS): beginning transition to session DAQ interface
%--------------------------------------------------------------------------
disp('nidaq_ai_init: starting NI hardware...');

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Initializes the NI board 
%------------------------------------------------------------------------
%------------------------------------------------------------------------


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Legacy Interface
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if strcmpi(iface, 'NI')
	%------------------------------------------------------------------------
	% 'nidaq' specifies the national instruments device with traditional
	% DAQ Toolbox interface, Device number 1 (get this from the 
	% NI Measurement & Automation Explorer (a.k.a., MAX) program)
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% CONFIGURE ANALOG INPUT SUBSYSTEM
	%------------------------------------------------------------------------
	fprintf('Initializing NIDAQ device for analog input...')
	try
		ai = analoginput('nidaq', Dnum);
		fprintf('...done\n')
	catch errEvent
		fprintf('\nProblem while initializing NIDAQ device!\n\n')
		disp(errEvent)
		return
	end
	% create AI channel
	fprintf('creating analog input channel...')
	chI = addchannel(ai, [0 1]);
	fprintf('...done\n');
	ai.Channel(1).ChannelName = 'responseL';
	ai.Channel(2).ChannelName = 'responseR';
	%-------------------------------------------------------
	% save in NI struct
	%-------------------------------------------------------
	NI.ai = ai;
	NI.chI = chI;
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
	% save in NI struct
	%---------------------------------------------------------------------
else
	errordlg('Unknown interface')
	error('Unknown interface');
end

