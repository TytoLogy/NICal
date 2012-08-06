%--------------------------------------------------------------------------
% NICal_caldata_init.m
%--------------------------------------------------------------------------
%	Script for NICal program to initialize/allocate caldata
%	structure for speaker calibration
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 1 March, 2012, branched from HeadphoneCal
%
% Revisions:
%	28 Mar 2012 (SJS):
%		-	made modifications in case FRdata are unused (i.e., calibration
%			mic is used instead of calibrated mic)
%	9 Jul 2012 (SJS) renamed for NICal project
%	12 Jul 2012 (SJS)
%	 -	added comments
%	 -	added leak variables for measuring crosstalk
%	2 Aug 2012 (SJS)
%	 - added background variable
%--------------------------------------------------------------------------

%------------------------------------------------------------
%------------------------------------------------------------
% Setup data storage variables and paths
%------------------------------------------------------------
%------------------------------------------------------------
caldata.time_str = datestr(now, 31);			% date and time
caldata.timestamp = now;							% timestamp
caldata.adFc = iodev.Fs;							% analog input rate
caldata.daFc = iodev.Fs;							% analog output rate
caldata.nrasters = Nfreqs;							% number of freqs to collect
caldata.range = [min(Freqs) max(Freqs)];		% freq range (matlab string)
caldata.reps = cal.Nreps;							% reps per frequency
caldata.settings = cal;								% cal struct
caldata.atten = cal.StartAtten;					% initial attenuator setting
caldata.max_spl = cal.Maxlevel;					% maximum spl
caldata.min_spl = cal.Minlevel;					% minimum spl

%------------------------------------------------------------------
%------------------------------------------------------------------
% initialize the caldata structure arrays for the calibration data
%------------------------------------------------------------------
%------------------------------------------------------------------
tmparr = zeros(handles.Nchannels, Nfreqs);
caldata.freq = Freqs;
caldata.mag = tmparr;
caldata.phase = tmparr;
caldata.dist = tmparr;
caldata.mag_stderr = tmparr;
caldata.phase_stderr = tmparr;
caldata.background = tmparr;
caldata.background_stderr = tmparr;

% if leak is to be measured, create storage space
if handles.cal.MeasureLeak
	caldata.leakmag = tmparr;
	caldata.leakmad_stderr = tmparr;
	caldata.leakphase = tmparr;
	caldata.leakphase_stderr = tmparr;
	caldata.leakdist = tmparr;
	caldata.leakdist_stderr = tmparr;
	caldata.leakdistphis = tmparr;
	caldata.leakdistphis_stderr = tmparr;
end

% if a non-reference mic (e.g., B&K) is not used, add some other stuff
if ~handles.cal.FRenable
	caldata.Side = handles.cal.Side;
	caldata.InputChannel = handles.cal.InputChannel;
	caldata.DAscale = handles.cal.DAscale;
	caldata.Gain_dB = handles.cal.MicGain;
	caldata.Gain = invdb(caldata.Gain_dB);
	caldata.CalMic_sense = handles.cal.MicSensitivity;
	caldata.VtoPa = VtoPa;
	caldata.frfile = '';
	caldata.frdata = [];
end
	
%------------------------------------------------------------------
%------------------------------------------------------------------
% Fetch the l and r headphone mic adjustment values for the 
% calibration frequencies using interpolation 
%------------------------------------------------------------------
%------------------------------------------------------------------
if ~handles.cal.FRenable
	frdata.lmagadjval = ones(size(caldata.freq));
	frdata.rmagadjval = ones(size(caldata.freq));
	frdata.lphiadjval = zeros(size(caldata.freq));
	frdata.rphiadjval = zeros(size(caldata.freq));
	frdata.DAscale = DAscale;
else
	frdata.lmagadjval = interp1(frdata.freq, frdata.ladjmag, caldata.freq);
	frdata.rmagadjval = interp1(frdata.freq, frdata.radjmag, caldata.freq);
	frdata.lphiadjval = interp1(frdata.freq, frdata.ladjphi, caldata.freq);
	frdata.rphiadjval = interp1(frdata.freq, frdata.radjphi, caldata.freq);
	caldata.frfile = handles.cal.mic_fr_file;
	caldata.DAscale = frdata.DAscale;
	caldata.frdata = frdata;
end


if DEBUG
	magsdbug = mags;
	phisdbug = phis;
end
	