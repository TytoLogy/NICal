%--------------------------------------------------------------------------
% HeaphoneCal_tdtinit.m
%--------------------------------------------------------------------------
% sets up TDT parameters
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sharad.shanbhag@einstein.yu.edu
%--------------------------------------------------------------------------
% Created: 5 Feb 2008
% 				Created from FFCal_tdtinit.m
% 
% Revisions:
%	19 June, 2009 (SJS): updated documentation
%	5 March, 2010 (SJS):	fixed status and tagname issue to comply
%								with new RPclose and RPgettag methods
%--------------------------------------------------------------------------

disp('...starting TDT hardware...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the TDT devices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Initialize RX8 device 2
	tmpdev = RX8init('GB', iodev.Dnum);
	iodev.C = tmpdev.C;
	iodev.handle = tmpdev.handle;
	iodev.status = tmpdev.status;
	
	% Initialize PA5 attenuators (left = 1 and right = 2)
	PA5L = PA5init('GB', 1);
	PA5R = PA5init('GB', 2);
	
	if cal.CheckCal
		iodev.REF = 3;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads circuits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	iodev.rploadstatus = RPload(iodev);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starts Circuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	RPrun(iodev);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	iodev.status = RPcheckstatus(iodev);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the tags and values for the circuit
% (added 5 Mar 2010 (SJS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tmptags = RPtagnames(iodev);
	iodev.TagName = tmptags;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Query the sample rate from the circuit and set up the time vector and 
% stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	iodev.Fs = RPsamplefreq(iodev);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up some of the buffer/stimulus parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% size of the Serial buffer
	npts=150000;  
	dt = 1/iodev.Fs;
	mclock=RPgettag(iodev, 'mClock');

	% Set the total sweep period time
	RPsettag(iodev, 'SwPeriod', ms2samples(cal.SweepPeriod, iodev.Fs));
	% Set the sweep count (may not be necessary)
	RPsettag(iodev, 'SwCount', 1);
	% Set the Stimulus Delay
	RPsettag(iodev, 'StimDelay', ms2samples(cal.StimDelay, iodev.Fs));
	% Set the Stimulus Duration
	RPsettag(iodev, 'StimDur', ms2samples(cal.StimDuration, iodev.Fs));
	% Set the length of time to acquire data
	RPsettag(iodev, 'AcqDur', ms2samples(cal.AcqDuration, iodev.Fs));
	% set the ttl pulse duration
	RPsettag(iodev, 'PulseDur', ms2samples(cal.TTLPulseDur, iodev.Fs));

	%Setup filtering if desired
	if cal.InputFilter
		RPsettag(iodev, 'HPFreq', cal.InputFc);
	end
	
	TDTINIT = 1;


