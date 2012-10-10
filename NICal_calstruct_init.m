function cal = NICal_calstruct_init()
%--------------------------------------------------------------------------
% cal = NICal_calstruct_init
%--------------------------------------------------------------------------
%	returns initialized cal structure that contains settings
% for running calibration in NICal
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 13 July, 2012
%
% Revisions:
% 	2 Aug 2012 (SJS): added CollectBackground
%	17 Aug 2012 (SJS): added FieldType for mic calibration, Nchannels
%	24 Aug 2012 (SJS): added TriggerTimeout, TriggerLevel for 
% 							 use in Triggered acquisition mode
%	27 Aug 2012 (SJS): modified fr and cal filenames
%	28 Aug 2012 (SJS): added gain, adjusted MicGain for Nchannels
%	10 Oct 2012 (SJS): added SaveRawData option
%--------------------------------------------------------------------------

%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
% Load Constants
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
NICal_Constants;

%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
% CALIBRATION SETTINGS
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------

%----------------------------------------------------------------
% set the 'side' to left channel (calibrate speaker on channel0)
%	from NICal_Constants:
%		left = 1
%		right = 2
%		both = 3
%----------------------------------------------------------------
cal.Side = L;

%----------------------------------------------------------------
% # reps per frequency
%----------------------------------------------------------------
cal.Nreps = 3;
%----------------------------------------------------------------

%----------------------------------------------------------------
%----------------------------------------------------------------
% Frequency range
%----------------------------------------------------------------
%----------------------------------------------------------------
cal.UseFreqList = 0;
% min, max, stepsize for frequencies (Hz)
cal.Fmin = 3000;
cal.Fmax = 3500;
cal.Fstep = 100;
% update Freqs and Nfreqs
cal.Freqs = cal.Fmin:cal.Fstep:cal.Fmax;
cal.Nfreqs = length(cal.Freqs);
%----------------------------------------------------------------

%----------------------------------------------------------------
%----------------------------------------------------------------
% DISPLAY/Plot
%----------------------------------------------------------------
%----------------------------------------------------------------
% decimation factor for plots
cal.deciFactor = 10;
%----------------------------------------------------------------

%----------------------------------------------------------------
%----------------------------------------------------------------
% ATTENUATION
%----------------------------------------------------------------
%----------------------------------------------------------------
% these are used to use a fixed attenuation level
% will essentially generate a frequency response curve
% for the speaker (with the microphone correction factor
% from *_fr.mat file from CalibrateHeadphoneMic applied)
%----------------------------------------------------------------
cal.AttenFix = 0;
cal.AttenFixValue = 90;
% set the min and max allowable stimulus levels (db SPL)
cal.Minlevel = 60;
cal.Maxlevel = 70;
% set the stepsize for adjusting attenuation (dB)
cal.AttenStep = 2;
% Set the starting attenuation value (dB)
% (better to set too high instead of too low!!!!)
cal.StartAtten = 90;

%----------------------------------------------------------------
% set the CheckCal flag in order to use a reference microphone to 
% check the calibration.  
%----------------------------------------------------------------
% 0 == no check, 1 = Left, 2 = Right, 3 = Both
%----------------------------------------------------------------
cal.CheckCal = 0;

%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
% INPUT/OUTPUT SETTINGS
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------

%----------------------------------------------------------------
% # channels
%----------------------------------------------------------------
cal.Nchannels = 2;

%---------------------------------------------
%---------------------------------------------
% set the stimulus/acquisition settings
%---------------------------------------------
%---------------------------------------------
% Stimulus Interval (ms)
cal.ISI = 100;
% Stimulus Duration (ms)
cal.StimDuration = 250;
% Delay of stimulus (ms)
cal.StimDelay = 5;
% Stimulus ramp on/off time (ms)
cal.StimRamp = 1;
% Duration of epoch (ms)
cal.SweepDuration = 500;
% some other derived variables (not user-settable)
% Total time to acquire data (ms)
%cal.AcqDuration = cal.SweepDuration;
% Total sweep time = sweep duration + inter stimulus interval (ms)
cal.SweepPeriod = cal.SweepDuration + cal.ISI;
%----------------------------------------------------------------
% output tone signal peak level
%----------------------------------------------------------------
cal.DAscale = 1;
%----------------------------------------------------------------
% Measure channel crosstalk (leak)
%----------------------------------------------------------------
cal.MeasureLeak = 0;
%----------------------------------------------------------------
% Save Raw Data
%----------------------------------------------------------------
cal.SaveRawData = 0;


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Define a bandpass filter for processing the data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%----------------------------------------------------------------
% Input Bandpass Filter (in software)
%----------------------------------------------------------------
cal.InputFilter = 1;
cal.InputHPFc = 100;
cal.InputLPFc = 125000;
% for NI card, set default sampling rate here (kludge!!!!)
cal.Fs = 500000;
% Nyquist frequency
fnyq = cal.Fs/2;
% filter order
cal.forder = 5;
% passband definition
cal.fband = [cal.InputHPFc cal.InputLPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[cal.fcoeffb, cal.fcoeffa] = butter(cal.forder, cal.fband, 'bandpass');
%----------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% triggered IO settings
%------------------------------------------------------------------------
% Triggered Acquisition bypasses most of the 
% settings and uses external TTL pulse to start
% sweep
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% enable (1)/disable (0) triggered Acquisition
cal.TriggeredAcquisition = 0;
%----------------------------------------------------------------
% *these settings only apply for triggered acquisition*
%----------------------------------------------------------------

%--------------------------
% trigger timeout
%--------------------------
% time to wait for TTL trigger in seconds
cal.TriggerSettings.TriggerTimeout = 10;

%--------------------------
% trigger level (TriggerConditionLevel)
%--------------------------
% triggering level (Volts) (used to set the TriggerConditionLevel setting
% in NI-DAQ hardware)
cal.TriggerSettings.TriggerLevel = 4;

%--------------------------
% trigger type
%--------------------------
% <string>, options:
% 'HwDigital'				The trigger source is an external digital signal. 
%	 							Pretrigger data cannot be captured. Control the trigger 
%								source with HwDigitalTriggerSource property. 
%								Specify the external digital signal with the 
%								TriggerCondition and TriggerConditionValue properties.
% 
% 'HwAnalogChannel'		The trigger source is an external analog signal 
%								(AI only). To set the trigger source, see 
%								TriggerChannel property.
% 
% 'HwAnalogPin'			The trigger source is a low-range external analog 
%								signal (AI only). Note that HwAnalogPin is supported 
%								only for Traditional NIDAQ devices. It is not supported 
%								for NIDAQmx devices.HWDigitalTriggerSource	
cal.TriggerSettings.TriggerType = 'HwDigital';

%--------------------------
% trigger source
%--------------------------
% if TriggerType = HwAnalogChannel or HwAnalogPin
% 	source for analog trigger input is first active	analog channel 
%	(integer from 0 to 15)
% else if TriggerType == HWDigitalTriggerSource: 
%	source for Digital trigger input
%	 'PFI0 to PFI15'			Use specified pin from PFI0 through PFI15.
%	 'RTSI0 to RTSI6'			Use specified pin from RTSI0 through RTSI6.
cal.TriggerSettings.TriggerSource = 'PFI0';

%--------------------------
% Trigger condition
%--------------------------
% The following trigger conditions are available for AI objects 
% when TriggerType is HwDigital. 
% 
% 'PositiveEdge'		The trigger occurs when the positive (rising) edge 
%								of a digital signal is detected.
% '{NegativeEdge}'	The trigger occurs when the negative (falling) edge 
%								of a digital signal is detected.
% 
% The following trigger conditions are available when TriggerType is 
% HwAnalogChannel or HwAnalogPin.
% 
% '{AboveHighLevel}'	The trigger occurs when the analog signal is above 
% 								the specified value.
% 'BelowLowLevel'		The trigger occurs when the analog signal is below 
% 								the specified value.
% 'InsideRegion'		The trigger occurs when the analog signal is inside 
% 								the specified region.
% 'LowHysteresis'		The trigger occurs when the analog signal is less 
% 								than the specified low value with hysteresis given 
%								by the specified high value.
% 'HighHysteresis'	The trigger occurs when the analog signal is 
% 								greater than the specified high value with 
% 								hysteresis given by the specified low value
cal.TriggerSettings.TriggerCondition = 'PositiveEdge';

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% MICROPHONE SETTINGS
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% use frequency response data when using a non-reference microphone
cal.FRenable = 0;
%----------------------------------------------------------------
% set the InputChannel to left channels;
%----------------------------------------------------------------
cal.InputChannel = L;
%----------------------------------------------------------------
% mic gain in dB
%----------------------------------------------------------------
cal.MicGain = 0;
% fix # gain values if Nchannels doesn't match # of gain values
if cal.Nchannels  ~= length(cal.MicGain)
	cal.MicGain = cal.MicGain(1) .* ones(1, cal.Nchannels);
end
% convert dB to linear scale
cal.Gain = invdb(cal.MicGain);
%----------------------------------------------------------------
% mic sensitivity in V/Pa
%----------------------------------------------------------------
cal.MicSensitivity = 0.1;
%----------------------------------------------------------------
% recording type (1 = freefield, 2 = pressure field (closed field)
%----------------------------------------------------------------
cal.FieldType = 1;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% FILE SETTINGS
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%----------------------------------------------------------------
% Auto save ear_cal.mat in experiment calibration data dir
%----------------------------------------------------------------
cal.AutoSave = 1;
%----------------------------------------------------------------
%----------------------------------------------------------------
% default files
%----------------------------------------------------------------
%----------------------------------------------------------------
datadir = ['C:\Users\' getenv('USERNAME')];
% default fr response file for Knowles mics
cal.mic_fr_file = fullfile(datadir, 'default.fr');
% default output file
cal.calfile = fullfile(datadir, 'ear.cal');

cal.CollectBackground = 0;



