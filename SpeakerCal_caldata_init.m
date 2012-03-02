%--------------------------------------------------------------------------
% HeadphoneCal_caldata_init.m
%--------------------------------------------------------------------------
%	Script for HeadphoneCal program to initialize/allocate caldata
%	structure for headphone speaker calibration
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbha@aecom.yu.edu
%--------------------------------------------------------------------------
% Created: 5 Feb 2008, branched from FFCal_settings.m
%
% Revisions:
%
%	23 January, 2009 (SJS):
%		-	renamed file from HeadphoneCal_caldata.m to 
%			HeadphoneCal_caldata_init.m to more
%			clearly indicate function of script.
%	19 June, 2009 (SJS): 
% 		-	added documentation, moved allocation of mags,
%			phis, etc arrays to _RunCalibration.m.
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
	caldata.frdata = frdata;
	caldata.atten = cal.StartAtten;					% initial attenuator setting
	caldata.max_spl = cal.Maxlevel;					% maximum spl
	caldata.min_spl = cal.Minlevel;					% minimum spl
	caldata.frfile = handles.cal.mic_fr_file;

	% set up the arrays to hold the data
	if cal.CheckCal
		if cal.CheckCal == L || cal.CheckCal == R
			Nchannels = 3;
		elseif cal.CheckCal == BOTH
			Nchannels = 4;
		end
	else
		Nchannels = 2;
	end
	
	%initialize the caldata structure arrays for the calibration data
	tmpcell = cell(Nchannels, Nfreqs);
	tmparr = zeros(Nchannels, Nfreqs);
	caldata.freq = Freqs;
	caldata.mag = tmparr;
	caldata.phase = tmparr;
	caldata.dist = tmparr;
	caldata.mag_stderr = tmparr;
	caldata.phase_stderr = tmparr;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch the l and r headphone mic adjustment values for the 
% calibration frequencies using interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	frdata.lmagadjval = interp1(frdata.freq, frdata.ladjmag, caldata.freq);
	frdata.rmagadjval = interp1(frdata.freq, frdata.radjmag, caldata.freq);
	frdata.lphiadjval = interp1(frdata.freq, frdata.ladjphi, caldata.freq);
	frdata.rphiadjval = interp1(frdata.freq, frdata.radjphi, caldata.freq);

	caldata.DAscale = frdata.DAscale;
	
	if DEBUG
		magsdbug = mags;
		phisdbug = phis;
	end
	