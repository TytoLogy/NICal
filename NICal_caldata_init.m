%--------------------------------------------------------------------------
% SpeakerCal_caldata_init.m
%--------------------------------------------------------------------------
%	Script for SpeakerCal program to initialize/allocate caldata
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
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup data storage variables and paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caldata.time_str = datestr(now, 31);			% date and time
caldata.timestamp = now;							% timestamp
caldata.adFc = iodev.Fs;							% analog input rate
caldata.daFc = iodev.Fs;							% analog output rate
caldata.nrasters = Nfreqs;							% number of freqs to collect
caldata.range = F;									% freq range (matlab string)
caldata.reps = cal.Nreps;							% reps per frequency
caldata.settings = cal;
caldata.atten = cal.StartAtten;					% initial attenuator setting
caldata.max_spl = cal.Maxlevel;					% maximum spl
caldata.min_spl = cal.Minlevel;					% minimum spl

% set up the arrays to hold the data
Nchannels = 2;

%initialize the caldata structure arrays for the calibration data
tmpcell = cell(Nchannels, Nfreqs);
tmparr = zeros(Nchannels, Nfreqs);
caldata.freq = Freqs;
caldata.mag = tmparr;
caldata.phase = tmparr;
caldata.dist = tmparr;
caldata.mag_stderr = tmparr;
caldata.phase_stderr = tmparr;

if ~handles.FRenable
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
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch the l and r headphone mic adjustment values for the 
% calibration frequencies using interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~handles.FRenable
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
	