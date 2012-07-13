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
%---------------------------------------------
% Global Constants
%---------------------------------------------
%---------------------------------------------
NICal_Constants;

%---------------------------------------------
%---------------------------------------------
% Load Microphone calibration data
%---------------------------------------------
%---------------------------------------------
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
	
%---------------------------------------------
%---------------------------------------------
% settings
%---------------------------------------------
%---------------------------------------------
% plot decimation factor
deciFactor = 1;
% read in the gain on the mic preamp
Gain_dB = handles.cal.MicGain;
% convert dB to linear scale
Gain = invdb(Gain_dB);
% this is the sensitivity of the calibration mic in V / Pa
% if FR file is used, get sens. from there, 
if read_ui_val(handles.FRenableCtrl) == 1
	CalMic_sense = frdata.calsettings.CalMic_sense;
else
	% otherwise, use cal information
	CalMic_sense = handles.cal.MicSensitivity;
end
% pre-compute the V -> Pa conversion factor
VtoPa = (CalMic_sense^-1);
% precompute the volts -> RMS conversion factor for sinusoids (0.7071)
RMSsin = 1/sqrt(2);

%---------------------------------------------
%---------------------------------------------
% set up the calibration frequency range
%---------------------------------------------
%---------------------------------------------
Freqs = handles.cal.Fmin:handles.cal.Fstep:handles.cal.Fmax;
F = [handles.cal.Fmin handles.cal.Fstep handles.cal.Fmax];
Nfreqs = length(Freqs);

%---------------------------------------------
%---------------------------------------------
% make local copy of iodev TDT control struct
%---------------------------------------------
%---------------------------------------------
iodev = handles.iodev;

