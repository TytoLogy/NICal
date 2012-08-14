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
%	17 Aug 2012 (SJS): added FieldType for mic calibration
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
% Triggered Acquisition bypasses most of the 
% settings and uses external TTL pulse to start
% sweep
%----------------------------------------------------------------
cal.TriggeredAcquisition = 0;
%----------------------------------------------------------------

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
% Attenuation
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

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Define a bandpass filter for processing the data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%----------------------------------------------------------------
% Input Bandpass Filter (in software)
%----------------------------------------------------------------
cal.InputFilter = 1;
cal.InputHPFc = 200;
cal.InputLPFc = 120000;
% for NI card, set default sampling rate here
cal.Fs = 250000;
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
cal.mic_fr_file = fullfile(datadir, 'default_fr.mat');
% default output file
cal.calfile = fullfile(datadir, 'ear_cal.mat');

cal.CollectBackground = 0;



