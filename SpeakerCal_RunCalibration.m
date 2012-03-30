%--------------------------------------------------------------------------
% SpeakerCal_RunCalibration.m
%--------------------------------------------------------------------------
% Runs the headphone speaker calibration using the in situ Knowles
% microphones as the the calibration mics, corrected using data from
% MicrophoneCal program (earphone fr data)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	March 2012 from HeadphoneCal_RunCalibration,	SJS
%
% Revisions:
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 1;
R = 2;
REF = 3;
BOTH = 3;
REFL = 3;
REFR = 4;
MAX_ATTEN = 120;

DEBUG = 0;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization Scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the COMPLETE flag to 0
COMPLETE = 0;

% Load the settings and constants 
SpeakerCal_settings;

% save the GUI handle information
guidata(hObject, handles);

% make a local copy of the cal settings structure
cal = handles.cal;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some checks and balances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if handles.FRenable
	% is frequency in range of the fr data for the headphones?
	% check low freq limit
	if F(1) < frdata.range(1)
		warning([mfilename ': requested LF calibration limit is out of FR file bounds']);
		return
	end
	% check high freq limit
	if F(3) > frdata.range(3)
		warning([mfilename ': requested HF calibration limit is out of FR file bounds']);
		return
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start TDT things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the TDT circuits
SpeakerCal_tdtinit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup caldata struct for storing the calibration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SpeakerCal_caldata_init;
% set the FRANGE output scale value (usually 5 V)
FRANGE = caldata.DAscale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate some arrays that are used locally
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = zeros(Nfreqs, cal.Nreps);
tmpcell = cell(Nchannels, 1);
for i=1:Nchannels
	tmpcell{i} = tmp;
end
mags = tmpcell;
magsraw = tmpcell;
phis = tmpcell;
phisraw = tmpcell;
dists = tmpcell;
distphis = tmpcell;
atten = tmpcell;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the start and end bins for the calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_bin = ms2bin(cal.StimDelay + cal.StimRamp, iodev.Fs);
if ~start_bin
	start_bin = 1;
end
end_bin = start_bin + ms2bin(cal.StimDuration-cal.StimRamp, iodev.Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create null stimulus and time vector for plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zerostim = syn_null(cal.StimDuration, iodev.Fs, 0);
zerostim = downsample(zerostim, deciFactor);
dt = deciFactor/iodev.Fs;
outpts = length(zerostim);		
% time vector for plots
tvec = 1000*dt*(0:(outpts-1));
zerostim = [0 0];		
acqpts = ms2bin(cal.AcqDuration, iodev.Fs);
stim_start = ms2bin(cal.StimDelay, iodev.Fs);
stim_end = stim_start + outpts - 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup attenuation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cal.AttenFix && between(cal.AttenFixValue, 0, 120)
	Latten = cal.AttenFixValue;
	Ratten = cal.AttenFixValue;
else
	% set the adjustable starting attenuator values	
	Latten = cal.StartAtten;
	Ratten = cal.StartAtten;
	if ~between(cal.AttenFixValue, 0, 120)
		warning([mfilename ': AttenFixValue out of range, using default StartAtten value'])
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now initiate sweeps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stopFlag = 0;
rep = 1;
freq_index = 1;
%*******************************LOOP through the frequencies
for freq = F(1):F(2):F(3)
	% update the frequency display value
	update_ui_str(handles.FreqVal, sprintf('%d', freq));

	% if we're collecting check data, print the frequency on the
	% command line
	disp(['FREQ: ' num2str(freq) '...']);

	% if cal.Side is 1 or 3 (LEFT or BOTH), calibrate L channel
	% 		cal.Side == 1 is LEFT
	% 		cal.Side == 2 is RIGHT
	% 		cal.Side == 3 is BOTH channels, 
	if cal.Side == 1 || cal.Side == 3

		% synthesize the L sine wave;
		[S, stimspec.RMS, stimspec.phi] = syn_calibrationtone2(cal.StimDuration, iodev.Fs, freq, 0, 'L');
		% scale the sound
		S = FRANGE * S;
		% apply the sin^2 amplitude envelope to the stimulus
		S = sin2array(S, cal.StimRamp, iodev.Fs);
		% plot the array
		axes(handles.Lstimplot);	plot(tvec, downsample(S(1, :), deciFactor), 'g');
		axes(handles.Rstimplot);	plot(zerostim, 'r');

		%loop while figuring out the L attenuator value.
		if cal.AttenFix
			% no need to test attenuation but, 
			% do need to set the attenuators
			PA5setatten(PA5L, Latten);
			PA5setatten(PA5R, MAX_ATTEN);
			update_ui_str(handles.LAttentext, Latten);
			update_ui_str(handles.RAttentext, MAX_ATTEN);
			% set retry to 0 to skip testing
			retry = 0;
		else
			retry = 1;
		end

		while retry
			% need to set the attenuators - since this is for the
			% L channel, set the R channel attenuator to MAX attenuation
			PA5setatten(PA5L, Latten);
			PA5setatten(PA5R, MAX_ATTEN);
			% update ui
			update_ui_str(handles.LAttentext, Latten);
			update_ui_str(handles.RAttentext, MAX_ATTEN);

			% play the sound;
			[resp, indx] = handles.iofunction(iodev, S, acqpts);
			% determine the magnitude and phase of the response
			[lmag, lphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, freq);
			% adjust for the gain of the preamp and apply correction
			% factors for RMS and microphone calibration

			lmag = RMSsin * lmag / (Gain*frdata.lmagadjval(freq_index));
			% compute dB SPL
			lmagdB = dbspl(VtoPa*lmag);
			update_ui_str(handles.LVal, sprintf('%.4f', lmag));
			update_ui_str(handles.LSPL, sprintf('%.4f', lmagdB));

			% check to see if the channel amplitude is in bounds
			if lmagdB > cal.Maxlevel
				Latten = Latten + cal.AttenStep;
				% if at limit, peg the attenuator value to max attenuation
				if Latten > MAX_ATTEN
					Latten = MAX_ATTEN;
					disp('Latten is maxed out!');
					retry = 0;
				end
			elseif lmagdB < cal.Minlevel
				Latten = Latten - cal.AttenStep;
				if Latten <= 0
					Latten = 0;
					disp('Latten at minimum level!');
					retry = 0;
				end
			else
				retry = 0;
			end

			% plot the response
			axes(handles.Lmicplot);
			plot(downsample(resp{L}(stim_start:stim_end), deciFactor), 'g');
			axes(handles.Rmicplot);
			plot(downsample(resp{R}(stim_start:stim_end), deciFactor), 'r');
		end

		pause(0.001*cal.ISI);

		% now, collect the data for frequency FREQ, LEFT headphone
		for rep = 1:cal.Nreps
			% play the sound;
			[resp, indx] = handles.iofunction(iodev, S, acqpts);

			% determine the magnitude and phase of the response
			[lmag, lphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, freq);
			[ldistmag, ldistphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, 2*freq);		

			% compute harmonic distortion measures before 
			% applying corrections for the knowles mic response
			dists{L}(freq_index, rep) = ldistmag / lmag;

			% adjust for the gain of the preamp and apply correction
			% factors for RMS and microphone calibration
			lmag_adjusted = RMSsin * lmag / (Gain*frdata.lmagadjval(freq_index));

			% Store the values in the cell arrays for later averaging
			% (we'll do the averages later in order to save time while
			%  running the calibration curves)
			% adjust for the gain of the preamp and convert to Pascals
			mags{L}(freq_index, rep) = VtoPa*(lmag_adjusted);
			phis{L}(freq_index, rep) = lphi - frdata.lphiadjval(freq_index);

			% store distortion and leak values
			distphis{L}(freq_index, rep) = ldistphi - frdata.lphiadjval(freq_index);

			% store the attenuator setting - will need this to compute
			% maximum attainable SPL at this frequency
			atten{L}(freq_index, rep) = Latten;

			% if we are collecting "check" data using a reference
			% microphone (i.e., B & K calibration mic), we have a few
			% more things to do...
			if cal.CheckCal == L
				[tmpmag, tmpphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, freq);
				mags{REF}(freq_index, rep) = VtoPa * RMSsin * tmpmag;
				phis{REF}(freq_index, rep) = tmpphi;
				[tmpdistmag, tmpdistphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, 2*freq);
				dists{REF}(freq_index, rep) = tmpdistmag;
				distphis{REF}(freq_index, rep) = tmpdistphi;
				disp(sprintf('ref mag: %f    L mag: %f', dbspl(mags{REF}(freq_index, rep)), dbspl(mags{L}(freq_index, rep)) ));
			elseif cal.CheckCal == BOTH
				[tmpmag, tmpphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, freq);
				mags{REFL}(freq_index, rep) = VtoPa * RMSsin * tmpmag;
				phis{REFL}(freq_index, rep) = tmpphi;
				[tmpdistmag, tmpdistphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, 2*freq);
				dists{REFL}(freq_index, rep) = tmpdistmag;
				distphis{REFL}(freq_index, rep) = tmpdistphi;
				disp(sprintf('refL mag: %f    L mag: %f', dbspl(mags{REFL}(freq_index, rep)), dbspl(mags{L}(freq_index, rep)) ));
			end

			% if DEBUG is set, save the raw magnitude and phase values
			if DEBUG
				magsdbug{L}(freq_index, rep) = lmag;
				phisdbug{L}(freq_index, rep) = lphi;

				if cal.CheckCal == L
					magsdbug{REF}(freq_index, rep) = tmpmag;
					phisdbug{REF}(freq_index, rep) = tmpphi;
				elseif cal.CheckCal == BOTH
					magsdbug{REFL}(freq_index, rep) = tmpmag;
					phisdbug{REFL}(freq_index, rep) = tmpphi;
				end
			end

			% plot the response
			axes(handles.Lmicplot);
			plot(downsample(resp{L}(stim_start:stim_end), deciFactor), 'g');
			axes(handles.Rmicplot);
			plot(downsample(resp{R}(stim_start:stim_end), deciFactor), 'r');

			update_ui_str(handles.LVal, sprintf('%.4f', lmag));
			update_ui_str(handles.LSPL, sprintf('%.4f', dbspl(mags{L}(freq_index, rep))));
			update_ui_str(handles.RVal, '---');
			update_ui_str(handles.RSPL, '---');

			pause(0.001*cal.ISI);
		end
	end

	pause(0.001*cal.ISI);

	if cal.Side == 2 || cal.Side == 3
		%RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
		% synthesize the R sine wave;
		[S, stimspec.RMS, stimspec.phi] = syn_calibrationtone2(cal.StimDuration, iodev.Fs, freq, 0, 'R');
		% scale the sound
		S = FRANGE * S;
		% apply the sin^2 amplitude envelope to the stimulus
		S = sin2array(S, cal.StimRamp, iodev.Fs);
		% plot the array
		axes(handles.Lstimplot);	plot(zerostim, 'g');
		axes(handles.Rstimplot);	plot(tvec, downsample(S(2, :), deciFactor), 'r');

		%loop while figuring out the R attenuator value.
		if cal.AttenFix
			% no need to test attenuation but, 
			% do need to set the attenuators
			PA5setatten(PA5L, MAX_ATTEN);
			PA5setatten(PA5R, Ratten);
			update_ui_str(handles.LAttentext, MAX_ATTEN);
			update_ui_str(handles.RAttentext, Ratten);
			% set retry to 0 to skip testing
			retry = 0;
		else
			retry = 1;
		end

		while retry
			% need to set the attenuators
			PA5setatten(PA5L, MAX_ATTEN);
			PA5setatten(PA5R, Ratten);
			update_ui_str(handles.LAttentext, MAX_ATTEN);
			update_ui_str(handles.RAttentext, Ratten);

			% play the sound;
			[resp, indx] = handles.iofunction(iodev, S, acqpts);
			% determine the magnitude and phase of the response
			[rmag, rphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, freq);
			% adjust for the gain of the preamp and apply correction
			% factors for RMS and microphone calibration
			rmag = RMSsin * rmag / (Gain*frdata.rmagadjval(freq_index));
			% compute dB SPL
			rmagdB = dbspl(VtoPa*rmag);
			update_ui_str(handles.RVal, sprintf('%.4f', rmag));
			update_ui_str(handles.RSPL, sprintf('%.4f', rmagdB));

			% check to see if the channel amplitude is in bounds
			if rmagdB > cal.Maxlevel
				Ratten = Ratten + cal.AttenStep;
				% if we're at the limit, peg the attenuator value to
				% max attenuation
				if Ratten > MAX_ATTEN
					Ratten = MAX_ATTEN;
					disp('Ratten is maxed out!');
					retry = 0;
				end
			elseif rmagdB < cal.Minlevel
				Ratten = Ratten - cal.AttenStep;
				if Ratten <= 0
					Ratten = 0;
					disp('Ratten at minimum level!');
					retry = 0;
				end
			else
				retry = 0;
			end

			% plot the response
			axes(handles.Lmicplot);
			plot(downsample(resp{L}(stim_start:stim_end), deciFactor), 'g');
			axes(handles.Rmicplot);
			plot(downsample(resp{R}(stim_start:stim_end), deciFactor), 'r');
		end

		pause(0.001*cal.ISI);

		% now, collect the data for frequency FREQ, RIGHT headphone
		for rep = 1:cal.Nreps
			% play the sound;
			[resp, indx] = handles.iofunction(iodev, S, acqpts);

			% determine the magnitude and phase of the response
			[rmag, rphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, freq);
			[rdistmag, rdistphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, 2*freq);				

			% compute distortion measures before applying corrections
			dists{R}(freq_index, rep) = rdistmag / rmag;

			% adjust for the gain of the preamp and apply correction
			% factors for RMS and microphone calibration
			rmag_adjusted = RMSsin * rmag / (Gain*frdata.rmagadjval(freq_index));

			% convert to Pascals (rms) and adjust phase measurements
			mags{R}(freq_index, rep) = VtoPa*(rmag_adjusted);
			phis{R}(freq_index, rep) = rphi - frdata.rphiadjval(freq_index);
			distphis{R}(freq_index, rep) = rdistphi - frdata.rphiadjval(freq_index);

			% store attenuation
			atten{R}(freq_index, rep) = Ratten;

			if cal.CheckCal == R
				[tmpmag, tmpphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, freq);
				mags{REF}(freq_index, rep) = VtoPa * RMSsin * tmpmag;
				phis{REF}(freq_index, rep) = tmpphi;
				[tmpdistmag, tmpdistphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, 2*freq);
				dists{REF}(freq_index, rep) = tmpdistmag;
				distphis{REF}(freq_index, rep) = tmpdistphi;					
				disp(sprintf('ref mag: %f    R mag: %f', dbspl(mags{REF}(freq_index, rep)), dbspl(mags{R}(freq_index, rep)) ));
			elseif cal.CheckCal == BOTH
				[tmpmag, tmpphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, freq);
				mags{REFR}(freq_index, rep) = VtoPa * RMSsin * tmpmag;
				phis{REFR}(freq_index, rep) = tmpphi;
				[tmpdistmag, tmpphi] = fitsinvec(resp{handles.cal.InputChannel}(start_bin:end_bin), 1, iodev.Fs, 2*freq);
				dists{REFR}(freq_index, rep) = tmpdistmag;
				distphis{REFR}(freq_index, rep) = tmpdistphi;
				disp(sprintf('refR mag: %f    R mag: %f', dbspl(mags{REFR}(freq_index, rep)), dbspl(mags{R}(freq_index, rep)) ));
			end

			% if DEBUG is set, save the raw magnitude and phase values
			if DEBUG
				magsdbug{R}(freq_index, rep) = rmag;
				phisdbug{R}(freq_index, rep) = rphi;

				if cal.CheckCal == R
					magsdbug{REF}(freq_index, rep) = tmpmag;
					phisdbug{REF}(freq_index, rep) = tmpphi;
				elseif cal.CheckCal == BOTH
					magsdbug{REFR}(freq_index, rep) = tmpmag;
					phisdbug{REFR}(freq_index, rep) = tmpphi;
				end
			end

			% plot the response
			axes(handles.Lmicplot);
			plot(downsample(resp{L}(stim_start:stim_end), deciFactor), 'g');
			axes(handles.Rmicplot);
			plot(downsample(resp{R}(stim_start:stim_end), deciFactor), 'r');
			% update values in text fields
			update_ui_str(handles.LVal, '---');
			update_ui_str(handles.LSPL, '---');
			update_ui_str(handles.RVal, sprintf('%.4f', rmag));
			update_ui_str(handles.RSPL, sprintf('%.4f', dbspl(mags{R}(freq_index, rep))));

			pause(0.001*cal.ISI);
		end
	end

	if read_ui_val(handles.Abort_ctrl) == 1
		disp('abortion detected')
		break
	end

	freq_index = freq_index + 1;
end %********************End of Cal loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exit gracefully (close TDT objects, etc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SpeakerCal_tdtexit;

if freq == F(3)
	COMPLETE = 1;
else
	return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Averages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_index = 1;
%*******************************LOOP through the frequencies
for freq = F(1):F(2):F(3)
	% compute the averages for this frequency

	% magnitude (dB) = db(rms) + atten
	magsraw{L}(freq_index, :) = dbspl(mags{L}(freq_index, :));
	magsraw{R}(freq_index, :) = dbspl(mags{R}(freq_index, :));
	mags{L}(freq_index, :) = dbspl(mags{L}(freq_index, :)) + atten{L}(freq_index, :);
	mags{R}(freq_index, :) = dbspl(mags{R}(freq_index, :)) + atten{R}(freq_index, :);

	% if Check data, save it
	if cal.CheckCal == L
		magsraw{REF}(freq_index, :) = dbspl(mags{REF}(freq_index, :));
		mags{REF}(freq_index, :) = dbspl(mags{REF}(freq_index, :)) + atten{L}(freq_index, :);
	elseif cal.CheckCal == R
		magsraw{REF}(freq_index, :) = dbspl(mags{REF}(freq_index, :));
		mags{REF}(freq_index, :) = dbspl(mags{REF}(freq_index, :)) + atten{R}(freq_index, :);
	elseif cal.CheckCal == BOTH
		magsraw{REFL}(freq_index, :) = dbspl(mags{REFL}(freq_index, :));
		magsraw{REFR}(freq_index, :) = dbspl(mags{REFR}(freq_index, :));
		mags{REFL}(freq_index, :) = dbspl(mags{REFL}(freq_index, :)) + atten{L}(freq_index, :);
		mags{REFR}(freq_index, :) = dbspl(mags{REFR}(freq_index, :)) + atten{R}(freq_index, :);
	end

	for channel = 1:Nchannels				
		caldata.mag(channel, freq_index) = mean( mags{channel}(freq_index, :) );
		caldata.mag_stderr(channel, freq_index) = std( mags{channel}(freq_index, :) );

		caldata.phase(channel, freq_index) = mean( unwrap(phis{channel}(freq_index, :)) );
		caldata.phase_stderr(channel, freq_index) = std( unwrap(phis{channel}(freq_index, :)) );

		caldata.dist(channel, freq_index) = mean( dists{channel}(freq_index, :) );
		caldata.dist_stderr(channel, freq_index) = std( dists{channel}(freq_index, :) );
	end
	freq_index = freq_index + 1;
end

caldata.magsraw = magsraw;
caldata.atten = atten;

if DEBUG
	caldata.magsdbug = magsdbug;
	caldata.phisdbug = phisdbug;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save handles and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.caldata = caldata;

guidata(hObject, handles);

% plot the calibration data
PlotCal(caldata);

if cal.AutoSave
	disp(['Saving calibration data in ' earcalfile ' ...']);
	save(earcalfile, '-MAT', 'caldata');
end

disp('Finished.')


save test.mat resp indx