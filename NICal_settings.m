%--------------------------------------------------------------------------
% NICal_settings.m
%--------------------------------------------------------------------------
% This sets up the NICal parameters
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 1 March, 2012
% 				Created from HeadphoneCal_settings.m
%
% Revisions:
%	9 July 2012 (SJS)
% 	 -	modifications to use NI hardware
%--------------------------------------------------------------------------

disp('...general setup starting...');

%---------------------------------------------
% Global Constants
%---------------------------------------------
	L = 1;
	R = 2;
	REF = 3;
	BOTH = 3;
	MAX_ATTEN = 120;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Microphone calibration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if read_ui_val(handles.FRenableCtrl) == 0
		DAscale = read_ui_str(handles.DAscaleCtrl, 'n');
		handles.cal.mic_fr = [];
		cal.DAscale = DAscale;
	else
		load(handles.cal.mic_fr_file, 'frdata');
		if ~isfield(frdata, 'DAscale')
			frdata.DAscale = frdata.calsettings.DAscale;
		end
		handles.cal.mic_fr = frdata;
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	earcalpath = pwd;
	earcalfile = fullfile(earcalpath, 'ear_cal.mat');

	deciFactor = 1;

	% read in the gain on the mic preamp
	Gain_dB = handles.cal.MicGain;

	Gain = invdb(Gain_dB);
	
	% this is the sensitivity of the calibration mic in V / Pa
	if read_ui_val(handles.FRenableCtrl) == 1
		CalMic_sense = frdata.calsettings.CalMic_sense;
	else
		CalMic_sense = handles.cal.MicSensitivity;
	end
	
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
	handles.cal.InputHPFc = 120;
	handles.cal.InputLPFc = 45000;
	%TTL pulse duration (msec)
	handles.cal.TTLPulseDur = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make local copy of iodev TDT control struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	iodev = handles.iodev;

