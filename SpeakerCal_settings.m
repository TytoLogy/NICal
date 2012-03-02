%--------------------------------------------------------------------------
% HeadphoneCal_settings.m
%--------------------------------------------------------------------------
% This sets up the HeaphoneCal parameters
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbha@aecom.yu.edu
%--------------------------------------------------------------------------
% Revisions:
%
%	5 Feb 2008:	Created from FFCal_settings.m
%	23 January, 2009 (SJS):
%		-	renamed miclrcal, replaced with frdata struct
% 	19 June, 2009 (SJS): added documentation
%	2 Nov, 2010 (SJS): moved iodev (TDT device) configuration to 
% 						HeadphoneCal_OpeningFcn() for initialization and
% 						to TDTSettingsMenuCtrl() for user modification
%--------------------------------------------------------------------------

disp('...general setup starting...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	L = 1;
	R = 2;
	REF = 3;
	BOTH = 3;
	MAX_ATTEN = 120;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Microphone calibration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	load(handles.cal.mic_fr_file, 'frdata');
	if ~isfield(frdata, 'DAscale')
		frdata.DAscale = frdata.calsettings.DAscale;
	end
	handles.cal.mic_fr = frdata;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set global settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	earcalpath = pwd;
	earcalfile = [earcalpath '\ear_cal.mat'];

	deciFactor = 1;

	% read in the gain on the mic preamp
	Gain_dB = [40 40];
	Gain = 10.^(Gain_dB./20);

	% this is the sensitivity of the calibration mic in V / Pa
	CalMic_sense = frdata.calsettings.CalMic_sense;
	
	% pre-compute the V -> Pa conversion factor
	VtoPa = (CalMic_sense^-1);

	% precompute the volts -> RMS conversion factor for sinusoids (0.7071)
	RMSsin = 1/sqrt(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the stimulus/acquisition settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% set up the calibration frequency range
	Freqs = handles.cal.Fmin:handles.cal.Fstep:handles.cal.Fmax;
	F = [handles.cal.Fmin handles.cal.Fstep handles.cal.Fmax];
	Nfreqs = length(Freqs);

	% Stimulus Interval (ms)
	handles.cal.StimInterval = 0;
	% Stimulus Duration (ms)
	handles.cal.StimDuration = 100;
	% Duration of epoch (ms)
	handles.cal.SweepDuration = 120;
	% Delay of stimulus (ms)
	handles.cal.StimDelay = 5;
	% Total time to acquire data (ms)
	handles.cal.AcqDuration = handles.cal.SweepDuration;
	% Total sweep time = sweep duration + inter stimulus interval (ms)
	handles.cal.SweepPeriod = handles.cal.SweepDuration + handles.cal.StimInterval;
	% Stimulus ramp on/off time
	handles.cal.StimRamp = 5;

	%Input Filter Fc
	handles.cal.InputFilter = 1;
	handles.cal.InputFc = 120;
	%TTL pulse duration (msec)
	handles.cal.TTLPulseDur = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make local copy of iodev TDT control struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	iodev = handles.iodev;

