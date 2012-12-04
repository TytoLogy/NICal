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
% 	27 Aug 2012 (SJS): removed handles.deciFactor 
%							(changed to handles.cal.deciFactor in main function)
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
	handles.cal.DAscale = DAscale;
else
	load(handles.cal.mic_fr_file, 'frdata', '-MAT');
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
% fix # gain values if Nchannels doesn't match # of gain values
if handles.cal.Nchannels  ~= length(handles.cal.MicGain)
	handles.cal.MicGain = handles.cal.MicGain(1) .* ones(1, handles.cal.Nchannels);
	update_ui_str(handles.MicGainCtrl, handles.cal.MicGain);
	guidata(hObject, handles);
end
% convert dB to linear scale
Gain = invdb(handles.cal.MicGain);
% this is the sensitivity of the calibration mic in V / Pa
% if FR file is used, get sens. from there, 
if read_ui_val(handles.FRenableCtrl) == 1
	CalMic_sense = frdata.calsettings.MicSensitivity;
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
Freqs = handles.cal.Freqs;
Nfreqs = length(Freqs);

%-----------------------------------------------------------------------
% check calibration frequency range
%------------------------------------------------
if handles.cal.FRenable
	% is frequency in range of the fr data for the headphones?
	% check low freq limit
	if Freqs(1) < frdata.range(1)
		warning('NICal:Freq', [mfilename ': requested Min calibration frequency is out of FR file bounds']);
		return
	end
	% check high freq limit
	if Freqs(end) > frdata.range(2)
		warning('NICal:Freq', [mfilename ': requested Max calibration frequency is out of FR file bounds']);
		return
	end
end


