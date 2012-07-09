function NI = nidaq_init(iface, Dnum)
%--------------------------------------------------------------------------
% NI = nidaq_init.m
%--------------------------------------------------------------------------
% NICal program
% TytoLogy Project
% initializes nidaq system
%------------------------------------------------------------------------
% Input Arguments:
% 	iface		must be 'NI'
%	Dnum		device id (usually 'Dev1')
% 
% Output Arguments:
% 	NI		struct containing settings for requested type
% 		NI.ao		analog output object
% 		NI.ai		analog input object
% 		NI.chO	analog output channel object
% 		NI.chI	analog input channel object
%------------------------------------------------------------------------
% See also: NICal
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
%--------------------------------------------------------------------------

disp('...starting NI hardware...');

if ~strcmpi(iface, 'NI')
	error('%s: invalid interface %s', mfilename, iface);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Now, Initialize the NI board (PCIe-6351)
%------------------------------------------------------------------------
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

ai.Channel(1).ChannelName = 'output_signal';
ai.Channel(2).ChannelName = 'microphone_signal';

%------------------------------------------------------------------------
% CONFIGURE ANALOG OUTPUT SUBSYSTEM
%------------------------------------------------------------------------
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
fprintf('creating analog output channel...')
chO = addchannel(ao, [0 1]);
fprintf('...done\n');

ao.Channel(1).ChannelName = 'stimulusL';
ao.Channel(2).ChannelName = 'stimulusR';

%-------------------------------------------------------
% set logging mode
%	'Disk'	sets logging mode to a file on disk (specified by 'LogFileName)
%	'Memory'	sets logging mode to memory only
%	'Disk&Memory'	logs to file and memory
%-------------------------------------------------------
set(ai, 'LoggingMode', 'Memory');

%-------------------------------------------------------
% save in NI struct
%-------------------------------------------------------
NI.ai = ai;
NI.ao = ao;
NI.chI = chI;
NI.chO = chO;
NI.status = 1;



