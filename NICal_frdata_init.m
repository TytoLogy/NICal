%--------------------------------------------------------------------------
% NICal_frdata_init.m
%--------------------------------------------------------------------------
%	Script to initialize/allocate frdata structure that holds the microphone
%	frequency response calibration data for the earphone microphones
%--------------------------------------------------------------------------
% Data Format:
% 
% 	frdata fields: 
% 	
% 			FIELD       FMT		   DESCRIPTION	
%          time_str: (str)			date and time of data collection
%         timestamp: (dbl)			matlab timestamp at start of data acq.
%              adFc: (dbl)			analog-digital conversion rate for data
%              daFc: (dbl)			digital-analog conversion rate for signals
%          nrasters: (dbl)			# of frequencies tested
%             range: [1X3 dbl]		array of Fmin Fstep Fmax freqs
%              reps: (dbl)			# of reps at each frequency
%       calsettings: (struct)		calibration settings structure (see MicrophoneCal_settings.m)
%             atten: (dbl)			attenuator setting, dB
%           max_spl: (dbl)			maximum dB SPL level
%           min_spl: (dbl)			minimum dB SPL signal level
% 			 DAscale: (dbl)			scaling factor for output signal in Volts
%              freq: [1x473 dbl]	frequencies
%               mag: [3x473 dbl]	magnitude data, (Left Vrms, Right Vrms, Ref VRMS)
%             phase: [3x473 dbl]	phase data (degrees)
%              dist: [3x473 dbl]	mag. distortion data (2nd harmonic)
%        mag_stderr: [3x473 dbl]	magnitude std. error.
%      phase_stderr: [3x473 dbl]	phase std. error
%        background: [3x2 dbl]	Background RMS level, Volts, (L, R, Ref channels)
%     bkpressureadj: [1x473 dbl]	ref. mic correction factor for pressure field measurements
%           ladjmag: [1x473 dbl]	L channel microphone magnitude correction factor (Vrms/Pascal_rms)
%           radjmag: [1x473 dbl]	R channel microphone magnitude correction factor (Vrms/Pascal_rms)
%           ladjphi: [1x473 dbl]	L channel microphone phase correction factor (deg)
%           radjphi: [1x473 dbl]	R channel microphone phase correction factor (deg)
% 
% 		microphone adjust factors are:
% 			magnitude adj = knowles mic Vrms / Ref mic Vrms
% 			phase adj = Knowls mic deg - ref mic degrees
% 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Revisions:
%
% 16 August, 2012 (SJS):	Created from MicrophoneCal_frdata_init.m
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup data storage variables and paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	frdata.time_str = datestr(now, 31);		% date and time
	frdata.timestamp = now;						% timestamp
	frdata.adFc = handles.iodev.Fs;						% analog input rate
	frdata.daFc = handles.iodev.Fs;						% analog output rate
	frdata.nrasters = handles.cal.Nfreqs;					% number of freqs to collect
	frdata.range = [min(Freqs) max(Freqs)];								% freq range (matlab string)
	frdata.reps = handles.cal.Nreps;							% reps per frequency
	frdata.calsettings = handles.cal;					% parameters for calibration session
	frdata.atten = handles.cal.AttenFixValue;					% initial attenuator setting
	frdata.max_spl = 0;							% maximum spl (will be determined in program)
	frdata.min_spl = 0;							% minimum spl (will be determined in program)
	frdata.DAscale = handles.cal.DAscale;				% output peak voltage level

	% set up the arrays to hold the data
	tmp = zeros(handles.cal.Nfreqs, handles.cal.Nreps);
	tmpcell = cell(handles.cal.Nchannels, 1);
	background = tmpcell;
	for i=1:handles.cal.Nchannels
		tmpcell{i} = tmp;
		background{i} = zeros(1, handles.cal.Nreps);
	end
	mags = tmpcell;
	phis = tmpcell;
	dists = tmpcell;
	
	%initialize the frdata structure arrays for the calibration data
	tmpcell = cell(handles.cal.Nchannels, handles.cal.Nfreqs);
	tmparr = zeros(handles.cal.Nchannels, handles.cal.Nfreqs);
	frdata.freq = Freqs;
	frdata.mag = tmparr;
	frdata.phase = tmparr;
	frdata.dist = tmparr;
	frdata.mag_stderr = tmparr;
	frdata.phase_stderr = tmparr;
	frdata.background = zeros(handles.cal.Nchannels, 2);
	if handles.cal.FieldType == 2
		frdata.rawmags = tmp;
	end

	% setup cell for raw data 
	rawdata.background = cell(handles.cal.Nreps, 1);
	rawdata.resp = cell(handles.cal.Nfreqs, handles.cal.Nreps);
	rawdata.Freqs = Freqs;
