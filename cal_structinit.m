function cal = cal_structinit
%--------------------------------------------------------------------------
% cal = cal_structinit
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
%--------------------------------------------------------------------------

%----------------------------------------------------------------
% Frequency range
%----------------------------------------------------------------
cal.Fmin = 3000;
cal.Fstep = 100;
cal.Fmax = 3500;

%----------------------------------------------------------------
% Attenuation
%----------------------------------------------------------------
% set the min and max allowable stimulus levels
cal.Minlevel = 60;
cal.Maxlevel = 70;
% Set the starting attenuation value
% (better to set too high instead of too low!!!!)
cal.StartAtten = 90;
% set the stepsize for adjusting attenuation
cal.AttenStep = 2;
%----------------------------------------------------------------
% these are used to use a fixed attenuation level
% will essentially generate a frequency response curve
% for the speaker (with the microphone correction factor
% from *_fr.mat file from CalibrateHeadphoneMic applied)
%----------------------------------------------------------------
cal.AttenFix = 0;
cal.AttenFixValue = 90;

%----------------------------------------------------------------
% # reps per frequency
%----------------------------------------------------------------
cal.Nreps = 3;

%----------------------------------------------------------------
% Inter-stim interval in seconds
%----------------------------------------------------------------
cal.ISI = 200;

%----------------------------------------------------------------
% Auto save ear_cal.mat in experiment calibration data dir
%----------------------------------------------------------------
cal.AutoSave = 1;

%----------------------------------------------------------------
% set the 'side' to left channels;
%----------------------------------------------------------------
cal.Side = 1;
%----------------------------------------------------------------
% set the InputChannel to left channels;
%----------------------------------------------------------------
cal.InputChannel = 1;

%----------------------------------------------------------------
% set the CheckCal flag in order to use a reference microphone to 
% check the calibration.  
%----------------------------------------------------------------
% 0 == no check, 1 = Left, 2 = Right, 3 = Both
%----------------------------------------------------------------
cal.CheckCal = 0;

%----------------------------------------------------------------
%----------------------------------------------------------------
% default files
%----------------------------------------------------------------
%----------------------------------------------------------------
% default fr response file for Knowles mics
cal.mic_fr_file = '..\CalibrationData\FFamp_CIThp_24-Sep-2009_fr.mat';
% default output file
cal.calfile = fullfile(pwd, 'ear_cal.mat');


%----------------------------------------------------------------
% mic gain in dB
%----------------------------------------------------------------
cal.MicGain = 0;

%----------------------------------------------------------------
% mic sensitivity in V/Pa
%----------------------------------------------------------------
cal.MicSensitivity = 0.1;

%----------------------------------------------------------------
% output tone signal peak level
%----------------------------------------------------------------
cal.DAscale = 1;

%---------------------------------------------
%---------------------------------------------
% set the stimulus/acquisition settings
%---------------------------------------------
%---------------------------------------------
% Stimulus Interval (ms)
cal.StimInterval = 100;
% Stimulus Duration (ms)
cal.StimDuration = 250;
% Duration of epoch (ms)
cal.SweepDuration = 500;
% Delay of stimulus (ms)
cal.StimDelay = 5;
% Total time to acquire data (ms)
cal.AcqDuration = cal.SweepDuration;
% Total sweep time = sweep duration + inter stimulus interval (ms)
cal.SweepPeriod = cal.SweepDuration + cal.StimInterval;
% Stimulus ramp on/off time (ms)
cal.StimRamp = 1;

% for NI card, set default sampling rate here
cal.Fs = 250000;


% Input Filter Fc
cal.InputFilter = 1;
cal.InputHPFc = 200;
cal.InputLPFc = 120000;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Define a bandpass filter for processing the data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Nyquist frequency
fnyq = cal.Fs/2;
% filter order
cal.forder = 5;
% passband definition
cal.fband = [cal.InputHPFc cal.InputLPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[cal.fcoeffb, cal.fcoeffa] = butter(cal.forder, cal.fband, 'bandpass');

