%--------------------------------------------------------------------------
% NICal_RunCalibration.m
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% Runs the speaker calibration
%--------------------------------------------------------------------------
% if FR Correction is selected, apply mic correction using data from
% MicrophoneCal program (earphone fr data)
%--------------------------------------------------------------------------
% See also: NICal
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	March 2012 from HeadphoneCal_RunCalibration,	SJS
%
% Revisions:
%	9 July, 2012 (SJS) renamed for NICal project
%	2 July, 2012 (SJS): added background calculation
%	27 Aug 2012 (SJS): 
% 	 -	some comments added
%	 -	noticed that CheckCal is somewhat useless, so disabled it in GUI.
%	 -	should probably rework the way attenuation is done (manage both
%		Latten and Ratten instead of using MAX_ATTEN when needed), but 
%		this is low priority
%	 -	changed deciFactor to handles.cal.deciFactor
%	18 Jan 2017 (SJS): updated comments
%	1 Feb 2017 (SJS): updated for session interface
%	6 Feb 2017 (SJS): modifying stimulus to include pre/post stim time
%	9 Feb 2017 (SJS): some cleaning up, added check_output_file call
%--------------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Global Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_Constants;
% local settings
% set the COMPLETE flag to 0
COMPLETE = 0;
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Initialization Scripts
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%---------------------------------------------
% Load the settings and constants 
%---------------------------------------------
NICal_settings;
% save the GUI handle information
guidata(hObject, handles);
%-----------------------------------------------------------------------
% check output  file - if it exists, check with user
%-----------------------------------------------------------------------
calfile = check_output_file(handles);
if isequal(calfile, 0)
	return
else
	handles.cal.calfile = calfile;
	update_ui_str(handles.CalFileCtrl, handles.cal.calfile);
	guidata(hObject, handles);
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Start DAQ things
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
[handles, initFlag] = NICal_NIinit(handles);
guidata(hObject, handles);
if initFlag == 0
	warning('NICAL:HW', '%s: NIinit failure', mfilename)
	return
end
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Define a bandpass filter for processing the data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Nyquist frequency
fnyq = handles.iodev.Fs / 2;
% passband definition
handles.cal.fband = [handles.cal.InputHPFc handles.cal.InputLPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[handles.cal.fcoeffb, handles.cal.fcoeffa] = ...
					butter(handles.cal.forder, handles.cal.fband, 'bandpass');
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Setup raw data file if needed
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
if handles.cal.SaveRawData
	[pathstr, fname, fext] = fileparts(handles.cal.calfile);
	rawfile = fullfile(pathstr, [fname '.dat']);
	fp = fopen(rawfile, 'w');
	writeStruct(fp, handles.cal, 'cal');
	fclose(fp);
end
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Setup caldata struct for storing the calibration data
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_caldata_init;
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Preallocate some arrays that are used locally
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
tmp = zeros(Nfreqs, handles.cal.Nreps);
tmpcell = cell(handles.Nchannels, 1);
for i=1:handles.Nchannels
	tmpcell{i} = tmp;
end
mags = tmpcell;
magsraw = tmpcell;
phis = tmpcell;
phisraw = tmpcell;
dists = tmpcell;
distphis = tmpcell;
atten = tmpcell;
magsV = tmpcell;
% add leak information if needed
if handles.cal.MeasureLeak == 1
	leakmags = tmpcell;
	leakphis = tmpcell;
	leakdists = tmpcell;
	leakdistphis = tmpcell;
end
% if CollectBackground
if handles.cal.CollectBackground == 1
	bgmags = tmpcell;
	bgdb = bgmags;
end
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set the start and end bins for the calibration
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
[start_bin, end_bin] = startendbins(handles);
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% create null stimulus and time vector for plots, set up plots
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% calculate difference between SweepDuration and StimDelay + StimDuration
% this will be used to pad end of stimulus
% this is due to the session DAQ interface using # of cued output samples to
% determine the number of samples to read in. For the legacy interface,
% this shouldn't make a difference
PostDuration = handles.cal.SweepDuration - ...
						(handles.cal.StimDelay + handles.cal.StimDuration);
% make sure PostDuration is ok (greater than or equal to 0)
if PostDuration < 0
	errordlg('SweepDuration must be greater than StimDelay + StimDuration');
	NICal_NIexit;
	COMPLETE = 0;
	return
else
	% if ok, create poststim
	poststim = syn_null(PostDuration, handles.iodev.Fs, 0);
end
% create null stimulus
zerostim = syn_null(handles.cal.StimDuration, handles.iodev.Fs, 0);
% insert stim delay
zerostim = insert_delay(zerostim, handles.cal.StimDelay, handles.iodev.Fs);
% append post-stim
zerostim = [zerostim poststim];
% downsample (no need to plot all points)
zerostim = downsample(zerostim, handles.cal.deciFactor);
% downsample-factor adjusted sample interval
dt = handles.cal.deciFactor/handles.iodev.Fs;
% time vector for stimulus plots
tvec_stim = 1000*dt*(0:(length(zerostim)-1));
% fake acquired data
zeroacq = syn_null(handles.cal.SweepDuration, handles.iodev.Fs, 0);
zeroacq = downsample(zeroacq, handles.cal.deciFactor);
acqpts = length(zeroacq);
% time vector for stimulus plots
tvec_acq = 1000*dt*(0:(acqpts-1));
% compute # of points per sweep
SweepPoints = ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);
% stimulus start and end points
stim_start = ms2bin(handles.cal.StimDelay, handles.iodev.Fs);
stim_end = stim_start + ms2bin(handles.cal.StimDuration, handles.iodev.Fs) - 1;
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Build null output array
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% synthesize null stimulus;
Nullstim = syn_null(handles.cal.StimDuration, handles.iodev.Fs, 1);
% scale the sound
Nullstim = 0 * Nullstim;
% insert delay
Nullstim = insert_delay(Nullstim, handles.cal.StimDelay, handles.iodev.Fs);
% add end pad
Nullstim = [Nullstim syn_null(PostDuration, handles.iodev.Fs, 1)];
% downsampled, single channel version for plotting
Nullstim_downsample =  downsample(Nullstim(1, :), handles.cal.deciFactor);
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Setup Plots
%-----------------------------------------------------------------------
% to speed up plotting, the vectors Lacq, Racq, tvec_acq, L/Rfft, fvec
% are pre-allocated and then those arrys are used as XDataSource and
% YDataSource for the respective plots
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-------------------------------------------------------
% create arrays for plotting and plot them
%-------------------------------------------------------
% stim
Lstim = zerostim;
Rstim = zerostim;
% acq
Lacq = zeroacq;
Racq = zeroacq;
% FFT
nfft = length(start_bin:end_bin);
tmp = zeros(1, nfft);
[fvec, Lfft] = daqdbfft(tmp, handles.iodev.Fs, nfft);
Rfft = Lfft;
% convert fvec to kHz
fvec = 0.001 * fvec;
clear tmp
%-------------------------------------------------------
% plot null data, save handles for time-domain plots
%-------------------------------------------------------
% stimulus
H.Lstim = plot(handles.Lstimplot, tvec_stim, Lstim, 'g');
set(H.Lstim, 'XDataSource', 'tvec_stim', 'YDataSource', 'Lstim');
ylabel(handles.Lstimplot, 'V');
H.Rstim = plot(handles.Rstimplot, tvec_stim, Rstim, 'r');
set(H.Rstim, 'XDataSource', 'tvec_stim', 'YDataSource', 'Rstim');
% response
H.Lacq = plot(handles.Lmicplot, tvec_acq, Lacq, 'g');
set(H.Lacq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Lacq');
ylabel(handles.Lmicplot, 'V')
H.Racq = plot(handles.Rmicplot, tvec_acq, Racq, 'r');
set(H.Racq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Racq');
%-------------------------------------------------------
% plot null data, save handles for frequency-domain plots
%-------------------------------------------------------
H.Lfft = plot(handles.Lfftplot, fvec, Lfft);
set(H.Lfft, 'XDataSource', 'fvec', 'YDataSource', 'Lfft');
xlabel(handles.Lfftplot, 'Frequency (kHz)')
ylabel(handles.Lfftplot, 'dBV')
H.Rfft = plot(handles.Rfftplot, fvec, Rfft);
set(H.Rfft, 'XDataSource', 'fvec', 'YDataSource', 'Rfft');
xlabel(handles.Rfftplot, 'Frequency (kHz)');
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% setup attenuation
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% check if AttenFix is set and that it is in the range of [0, MAX_ATTEN]
if handles.cal.AttenFix && between(handles.cal.AttenFixValue, 0, MAX_ATTEN)
	Latten = handles.cal.AttenFixValue;
	Ratten = handles.cal.AttenFixValue;
else
	% use StartAtten
	if ~between(handles.cal.AttenFixValue, 0, MAX_ATTEN)
		warning('NICal:Atten', [mfilename ...
					': AttenFixValue out of range, using default StartAtten value'])
	end
	% set the adjustable starting attenuator values	
	Latten = handles.cal.StartAtten;
	Ratten = handles.cal.StartAtten;
end
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Now initiate sweeps
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
stopFlag = 0;
rep = 1;
freq_index = 1;
%---------------------------------------------------
% make a local copy of the cal settings structure
%---------------------------------------------------
cal = handles.cal;
%---------------------------------------------------
% make local copy of iodev  control struct
%---------------------------------------------------
iodev = handles.iodev;
%---------------------------------------------------
% measure!
%---------------------------------------------------
%*******************************LOOP through the frequencies
for F = 1:Nfreqs
	% assign current frequency from Freqs list to freq
	freq = Freqs(F);
	% update the frequency display value
	update_ui_str(handles.FreqValText, sprintf('%d', round(freq)));
	% check for abort button press
	if read_ui_val(handles.AbortCtrl) == 1
		% if so, stop
		disp('abortion detected')
		break
	end
	% if we're collecting check data, print the frequency on the
	% command line
	if cal.CheckCal
		disp(['FREQ: ' num2str(freq) '...']);
	end
	%------------------------------------------------------------------
	% if cal.Side is 1 or 3 (LEFT or BOTH), calibrate L channel
	% Side is set by the SideCtrl pulldown under CalibrationSettings
	% 		cal.Side == 1 is LEFT
	% 		cal.Side == 2 is RIGHT
	% 		cal.Side == 3 is BOTH channels, 
	%------------------------------------------------------------------
	if cal.Side == 1 || cal.Side == 3
		%LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
		% set input channel: if InputChannel is set to 3 (both), 
		% use left input when calibrating Left side
		%LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
		if handles.cal.InputChannel == 3
			inChan = 1;
		else
			inChan = handles.cal.InputChannel;
		end
		%-------------------------------------------------------
		% Build stimulus output array for this frequency
		%-------------------------------------------------------
		% synthesize the L sine wave;
		[S, stimspec.RMS, stimspec.phi] = ...
									syn_calibrationtone2(cal.StimDuration, ...
																	iodev.Fs, freq, 0, 'L');
		% scale the sound
		S = cal.DAscale * S;
		% apply the sin^2 amplitude envelope to the stimulus before adding 
		% pre and post zeros
		S = sin2array(S, cal.StimRamp, iodev.Fs);
		% insert delay, add zeros to pad end
		S = [insert_delay(S, cal.StimDelay, iodev.Fs) ...
									syn_null(PostDuration, iodev.Fs, 1)];
		% save in Satt
		Satt = S;
		% plot the stimuli - first downsample then set R stim to zero
		Lstim = downsample(S(1, :), handles.cal.deciFactor);
		Rstim = zerostim;
		refreshdata(H.Lstim, 'caller');
		refreshdata(H.Rstim, 'caller');
		%-------------------------------------------------------
		% loop while figuring out the L attenuator value.
		%-------------------------------------------------------
		if cal.AttenFix
			% no need to test attenuation but, 
			% do need to set the attenuators
			Satt(1, :) = handles.attfunction(S(1, :), Latten);
			Satt(2, :) = handles.attfunction(S(2, :), MAX_ATTEN);
			update_ui_str(handles.LAttenText, Latten);
			update_ui_str(handles.RAttenText, MAX_ATTEN);
			% set retry to 0 to skip testing
			retry = 0;
		else
			% set retry to 1 in order to enter the attenuation set loop
			retry = 1;
		end
		while retry
			% need to set the attenuators - since this is for the
			% L channel, set the R channel attenuator to MAX attenuation
			Satt(1, :) = handles.attfunction(S(1, :), Latten);
			Satt(2, :) = handles.attfunction(S(2, :), MAX_ATTEN);
			% update atten val display
			update_ui_str(handles.LAttenText, Latten);
			update_ui_str(handles.RAttenText, MAX_ATTEN);
			% play the sound;
			if handles.DAQSESSION
				[resp, indx] = handles.iofunction(iodev, Satt, ...
														handles.cal.SweepDuration);
			else
				[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
			end
			% filter the data if asked
			if handles.cal.InputFilter
				tmp = sin2array(resp{L}, 1, iodev.Fs);
				resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				tmp = sin2array(resp{R}, 1, iodev.Fs);
				resp{R} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				clear tmp
			end
			% determine the magnitude and phase of the response
			[lmag, lphi] = fitsinvec(resp{inChan}(start_bin:end_bin), ...
												1, iodev.Fs, freq);
			update_ui_str(handles.LValText, sprintf('%.4f', 1000*lmag));
			% adjust for the gain of the preamp and apply correction
			% factors for RMS and microphone calibration
			lmag = RMSsin * lmag / (Gain(L)*frdata.lmagadjval(freq_index));
			% compute dB SPL
			lmagdB = dbspl(VtoPa*lmag);
			update_ui_str(handles.LSPLText, sprintf('%.4f', lmagdB));
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
			Lacq = downsample(resp{L}, handles.cal.deciFactor);
			refreshdata(H.Lacq, 'caller');
			if handles.cal.MeasureLeak
				Racq = downsample(resp{R}, handles.cal.deciFactor);
				refreshdata(H.Racq, 'caller');
			end
		end		% END of L attenuation loop

		
		%-------------------------------------------------------
		% pause for ISI 
		%-------------------------------------------------------
		pause(0.001*cal.ISI);

		%-------------------------------------------------------
		% now, collect the data for frequency FREQ, LEFT channel
		%-------------------------------------------------------
		for rep = 1:cal.Nreps
			% update the reps display value
			update_ui_str(handles.RepNumText, sprintf('%d L', rep));
			% play the sound;
			if handles.DAQSESSION
				[resp, indx] = handles.iofunction(iodev, Satt, ...
														handles.cal.SweepDuration);
			else
				[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
			end
			% filter the data if asked
			if handles.cal.InputFilter
				tmp = sin2array(resp{L}, 1, iodev.Fs);
				resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				if handles.cal.MeasureLeak
					tmp = sin2array(resp{R}, 1, iodev.Fs);
					resp{R} = filtfilt(handles.cal.fcoeffb, ...
												handles.cal.fcoeffa, tmp);
				end
				clear tmp
			end

			% determine the magnitude and phase of the response
			[lmag, lphi] = fitsinvec(resp{inChan}(start_bin:end_bin), ...
											1, iodev.Fs, freq);
			[ldistmag, ldistphi] = ...
										fitsinvec(resp{inChan}(start_bin:end_bin), ...
															1, iodev.Fs, 2*freq);		
			update_ui_str(handles.LValText, sprintf('%.4f', 1000*lmag));

			% store peak mag (in Volts) in magsV 
			magsV{L}(freq_index, rep) = lmag;
			
			% compute harmonic distortion measures before 
			% applying corrections for the knowles mic response
			dists{L}(freq_index, rep) = ldistmag / lmag;

			% adjust for the gain of the preamp and apply correction
			% factors for RMS and microphone calibration
			lmag_adjusted = RMSsin * lmag / ...
									(Gain(L)*frdata.lmagadjval(freq_index));

			% Store the values in the cell arrays for later averaging
			% (we'll do the averages later in order to save time while
			%  running the calibration curves)
			% adjust for the gain of the preamp and convert to Pascals
			mags{L}(freq_index, rep) = VtoPa*(lmag_adjusted);
			phis{L}(freq_index, rep) = lphi - frdata.lphiadjval(freq_index);
			update_ui_str(handles.LSPLText, ...
								sprintf('%.4f', dbspl(mags{L}(freq_index, rep))));

			% store distortion and leak values
			distphis{L}(freq_index, rep) = ldistphi - ...
														frdata.lphiadjval(freq_index);

			% store the attenuator setting - will need this to compute
			% maximum attainable SPL at this frequency
			atten{L}(freq_index, rep) = Latten;

			% if we are collecting "check" data using a reference
			% microphone (i.e., B & K calibration mic), we have a few
			% more things to do...
			if cal.CheckCal == L
				% mag, phase of L channel
				[tmpmag, tmpphi] = ...
										fitsinvec(resp{inChan}(start_bin:end_bin), ...
															1, iodev.Fs, freq);
				% convert tmpmag to Pascals (rms)
				mags{REF}(freq_index, rep) = VtoPa * RMSsin * tmpmag;
				phis{REF}(freq_index, rep) = tmpphi;
				% distortion measures
				[tmpdistmag, tmpdistphi] = ...
											fitsinvec(resp{inChan}(start_bin:end_bin), ...
													1, iodev.Fs, 2*freq);
				dists{REF}(freq_index, rep) = tmpdistmag;
				distphis{REF}(freq_index, rep) = tmpdistphi;
				% feedback to user
				fprintf('ref mag: %f    L mag: %f', ...
										dbspl(mags{REF}(freq_index, rep)), ...
										dbspl(mags{L}(freq_index, rep)) );
			elseif cal.CheckCal == BOTH
				[tmpmag, tmpphi] = fitsinvec(resp{inChan}(start_bin:end_bin), ...
														1, iodev.Fs, freq);
				mags{REFL}(freq_index, rep) = VtoPa * RMSsin * tmpmag;
				phis{REFL}(freq_index, rep) = tmpphi;
				[tmpdistmag, tmpdistphi] = ...
										fitsinvec(resp{inChan}(start_bin:end_bin), 1, ...
														iodev.Fs, 2*freq);
				dists{REFL}(freq_index, rep) = tmpdistmag;
				distphis{REFL}(freq_index, rep) = tmpdistphi;
				fprintf('refL mag: %f    L mag: %f', ...
								dbspl(mags{REFL}(freq_index, rep)), ...
								dbspl(mags{L}(freq_index, rep)) );
			end

			% if DEBUG is set, save the raw magnitude and phase values
			if DEBUG
				magsdbug{L}(freq_index, rep) = lmag; %#ok<*SAGROW>
				phisdbug{L}(freq_index, rep) = lphi;

				if cal.CheckCal == L
					magsdbug{REF}(freq_index, rep) = tmpmag;
					phisdbug{REF}(freq_index, rep) = tmpphi;
				elseif cal.CheckCal == BOTH
					magsdbug{REFL}(freq_index, rep) = tmpmag;
					phisdbug{REFL}(freq_index, rep) = tmpphi;
				end
			end

			% if MeasureLeak is requested, measure it!
			if handles.cal.MeasureLeak
				% determine magnitude and phase of the response in the
				% opposite channel - this is the leak magnitude and phase
				[rleakmag, rleakphi] = ...
													fitsinvec(resp{R}(start_bin:end_bin), ...
																				1, iodev.Fs, freq);
				update_ui_str(handles.RValText, sprintf('%.4f', 1000*rleakmag));
				% compute leak distortion (1st harmonic)
				[rleakdistmag, rleakdistphi] = ...
													fitsinvec(resp{R}(start_bin:end_bin), ...
																			1, iodev.Fs, 2*freq);
				% compute harmonic distortion measures before 
				% applying corrections for the knowles mic response
				leakdists{R}(freq_index, rep) = rleakdistmag / rleakmag;
				% adjust for the gain of the preamp and apply correction
				% factors for RMS and microphone calibration
				rleakmag = RMSsin * rleakmag / ...
												(Gain(R)*frdata.rmagadjval(freq_index));
				% store leak values
				leakmags{R}(freq_index, rep) = VtoPa*(rleakmag);
				leakphis{R}(freq_index, rep) = ...
											rleakphi - frdata.rphiadjval(freq_index);
				leakdistphis{R}(freq_index, rep) = ...
											rleakdistphi - frdata.rphiadjval(freq_index);
				% update R text display
				update_ui_str(handles.RSPLText, sprintf('%.4f', ...
														dbspl(leakmags{R}(freq_index, rep))));
			else
				update_ui_str(handles.RValText, '---');
				update_ui_str(handles.RSPLText, '---');
			end
			
			% plot the response and FFT
			Lacq = downsample(resp{L}, handles.cal.deciFactor);
			refreshdata(H.Lacq, 'caller');
			[tmpf, Lfft] = daqdbfft(resp{L}(start_bin:end_bin), ...
											iodev.Fs, length(resp{L}(start_bin:end_bin)));
			refreshdata(H.Lfft, 'caller');
			if handles.cal.MeasureLeak
				Racq = downsample(resp{R}, handles.cal.deciFactor);
				refreshdata(H.Racq, 'caller');
				[tmpf, Rfft] = daqdbfft(resp{R}(start_bin:end_bin), ...
											iodev.Fs, length(resp{R}(start_bin:end_bin)));
				refreshdata(H.Rfft, 'caller');
			end
			drawnow
			% draw spectrogram
			axes(handles.Lspecgram); %#ok<*LAXES>
			myspectrogram(resp{L}, iodev.Fs, ...
									[10 5], @hamming, handles.SpectrumWindow, ...
									[-100 -1], false, 'default', false, 'per');
			if handles.cal.MeasureLeak
				axes(handles.Rspecgram);
				myspectrogram(resp{R}, iodev.Fs, ...
										[10 5], @hamming, handles.SpectrumWindow, ...
										[-100 -1], false, 'default', false, 'per');
			end
			% save raw data
			if handles.cal.SaveRawData
				fp = fopen(rawfile, 'a');
				if handles.cal.MeasureLeak
					writeCell(fp, resp);
				else
					writeCell(fp, resp(L));
				end
				fclose(fp);
			end
			% Pause for ISI
			pause(0.001*cal.ISI);
			% check for abort button press
			if read_ui_val(handles.AbortCtrl) == 1
				% if so, stop
				disp('abortion detected')
				break
			end
		end
		%---------------------------------------------------------------------
		% now, collect the background data for frequency FREQ, LEFT channel
		%---------------------------------------------------------------------
		if handles.cal.CollectBackground
			Lstim = Nullstim_downsample;
			Rstim = Nullstim_downsample;
			refreshdata(H.Lstim, 'caller');
			refreshdata(H.Rstim, 'caller');
			for rep = 1:cal.Nreps
				% update the reps display value
				update_ui_str(handles.RepNumText, sprintf('%d L (bg)', rep));
				% play the sound;
					if handles.DAQSESSION
						[resp, indx] = handles.iofunction(iodev, Satt, ...
																handles.cal.SweepDuration);
					else
						[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
					end
				% filter the data if asked
				if handles.cal.InputFilter
					tmp = sin2array(resp{L}, 1, iodev.Fs);
					resp{L} = filtfilt(handles.cal.fcoeffb, ...
															handles.cal.fcoeffa, tmp);
					tmp = sin2array(resp{R}, 1, iodev.Fs);
					resp{R} = filtfilt(handles.cal.fcoeffb, ...
															handles.cal.fcoeffa, tmp);
					clear tmp
				end
				% determine the magnitude and phase of the response
				bgmag = fitsinvec(resp{inChan}(start_bin:end_bin), ...
																			1, iodev.Fs, freq);
				% update LValText
				update_ui_str(handles.LValText, sprintf('%.4f', 1000*bgmag));
				% adjust for the gain of the preamp and apply correction
				% factors for RMS and microphone calibration
				bgmag_adjusted = RMSsin * bgmag / ...
										(Gain(L)*frdata.lmagadjval(freq_index));
				% Store the values in the cell arrays for later averaging
				% (we'll do the averages later in order to save time while
				%  running the calibration curves)
				% adjust for the gain of the preamp and convert to Pascals
				bgmags{L}(freq_index, rep) = VtoPa*(bgmag_adjusted);
				update_ui_str(handles.LSPLText, ...
										sprintf('%.4f', ...
														dbspl(bgmags{L}(freq_index, rep))));
				update_ui_str(handles.RValText, '---');
				update_ui_str(handles.RSPLText, '---');
				% plot the response
				Lacq = downsample(resp{L}, handles.cal.deciFactor);
				Racq = downsample(resp{R}, handles.cal.deciFactor);
				refreshdata(H.Lacq, 'caller');
				refreshdata(H.Racq, 'caller');
				% plot fft
				figure(10)
				[tmpf, tmpm] = daqdbfft(resp{L}(start_bin:end_bin), ...
										iodev.Fs, length(resp{L}(start_bin:end_bin)));
				plot(tmpf, tmpm);
				title('Left Background')
				% Pause for ISI
				if handles.cal.SaveRawData
					fp = fopen(rawfile, 'a');
					writeCell(fp, resp); 				
					fclose(fp);
				end
				pause(0.001*cal.ISI);
				% check for abort button press
				if read_ui_val(handles.AbortCtrl) == 1
					% if so, stop
					disp('abortion detected')
					break
				end
			end		
		end
	end		% END OF L CHANNEL	
	
	%------------------------------------
	% check for abort button press
	%------------------------------------
	if read_ui_val(handles.AbortCtrl) == 1
		% if so, stop
		disp('abortion detected')
		break
	end
	%------------------------------------
	% pause for ISI
	%------------------------------------
	pause(0.001*cal.ISI);

	%------------------------------------------------------------------
	% if cal.Side is 2 or 3 (RIGHT or BOTH), calibrate *R* channel
	%------------------------------------------------------------------	
	if cal.Side == 2 || cal.Side == 3
		%RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
		% set input channel: if INputChannel is set to 3 (both), use right
		if handles.cal.InputChannel == 3
			inChan = 2;
		else
			inChan = handles.cal.InputChannel;
		end
		% synthesize the R sine wave;
		[S, stimspec.RMS, stimspec.phi] = ...
					syn_calibrationtone2(cal.StimDuration, iodev.Fs, ...
													freq, 0, 'R');
		% scale the sound
		S = cal.DAscale * S;
		% apply the sin^2 amplitude envelope to the stimulus
		S = sin2array(S, cal.StimRamp, iodev.Fs);
		% insert delay, add zeros to pad end
		S = [insert_delay(S, cal.StimDelay, iodev.Fs) ...
									syn_null(PostDuration, iodev.Fs, 1)];
		% save in Satt
		Satt = S;
		% plot the stimulus arrays
		Lstim = zerostim;
		Rstim = downsample(S(2, :), handles.cal.deciFactor);
		refreshdata(H.Lstim, 'caller');
		refreshdata(H.Rstim, 'caller');
		% attenuation
		if cal.AttenFix
			% no need to test attenuation but, 
			% do need to set the attenuators
			Satt(1, :) = handles.attfunction(S(1, :), MAX_ATTEN);
			Satt(2, :) = handles.attfunction(S(2, :), Ratten);
			update_ui_str(handles.LAttenText, MAX_ATTEN);
			update_ui_str(handles.RAttenText, Ratten);
			% set retry to 0 to skip testing
			retry = 0;
		else
			retry = 1;
		end
		%loop while figuring out the R attenuator value.
		while retry
			% need to set the attenuators
			Satt(1, :) = handles.attfunction(S(1, :), MAX_ATTEN);
			Satt(2, :) = handles.attfunction(S(2, :), Ratten);
			update_ui_str(handles.LAttenText, MAX_ATTEN);
			update_ui_str(handles.RAttenText, Ratten);
			% play the sound;
			if handles.DAQSESSION
				[resp, indx] = handles.iofunction(iodev, Satt, ...
														handles.cal.SweepDuration);
			else
				[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
			end
			% filter the data if asked
			if handles.cal.InputFilter
				tmp = sin2array(resp{L}, 1, iodev.Fs);
				resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				tmp = sin2array(resp{R}, 1, iodev.Fs);
				resp{R} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				clear tmp
			end
			% determine the magnitude and phase of the response
			[rmag, rphi] = fitsinvec(resp{inChan}(start_bin:end_bin), 1, ...
												iodev.Fs, freq);
			% update text display
			update_ui_str(handles.RValText, sprintf('%.4f', 1000*rmag));
			% adjust for the gain of the preamp and apply correction
			% factors for RMS and microphone calibration
			rmag = RMSsin * rmag / (Gain(R)*frdata.rmagadjval(freq_index));
			% compute dB SPL
			rmagdB = dbspl(VtoPa*rmag);
			update_ui_str(handles.RSPLText, sprintf('%.4f', rmagdB));
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
			if handles.cal.MeasureLeak
				Lacq = downsample(resp{L}, handles.cal.deciFactor);
				refreshdata(H.Lacq, 'caller');
			end
			Racq = downsample(resp{R}, handles.cal.deciFactor);
			refreshdata(H.Racq, 'caller');
		end
		pause(0.001*cal.ISI);
		% now, collect the data for frequency FREQ, RIGHT headphone
		for rep = 1:cal.Nreps
			% update the reps display value
			update_ui_str(handles.RepNumText, sprintf('%d R', rep));
			% play the sound;
			if handles.DAQSESSION
				[resp, indx] = handles.iofunction(iodev, Satt, ...
														handles.cal.SweepDuration);
			else
				[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
			end
			% filter the data if asked
			if handles.cal.InputFilter
				if handles.cal.MeasureLeak
					tmp = sin2array(resp{L}, 1, iodev.Fs);
					resp{L} = filtfilt(handles.cal.fcoeffb, ...
													handles.cal.fcoeffa, tmp);
				end
				tmp = sin2array(resp{R}, 1, iodev.Fs);
				resp{R} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
				clear tmp
			end
			% determine the magnitude and phase of the response
			[rmag, rphi] = fitsinvec(resp{inChan}(start_bin:end_bin), ...
														1, iodev.Fs, freq);
			[rdistmag, rdistphi] = fitsinvec(resp{inChan}(start_bin:end_bin), ...
															1, iodev.Fs, 2*freq);				
			% store peak mag (in Volts) in magsV 
			magsV{R}(freq_index, rep) = rmag;
			% update values in text fields
			update_ui_str(handles.RValText, sprintf('%.4f', 1000*rmag));
			% compute distortion measures before applying corrections
			dists{R}(freq_index, rep) = rdistmag / rmag;
			% adjust for the gain of the preamp and apply correction
			% factors for RMS and microphone calibration
			rmag_adjusted = RMSsin * rmag / ...
										(Gain(R)*frdata.rmagadjval(freq_index));
			% update text display
			update_ui_str(handles.RSPLText, sprintf('%.4f', ...
																dbspl(VtoPa*rmag_adjusted)));
			% convert to Pascals (rms) and adjust phase measurements
			mags{R}(freq_index, rep) = VtoPa*(rmag_adjusted);
			phis{R}(freq_index, rep) = rphi - frdata.rphiadjval(freq_index);
			distphis{R}(freq_index, rep) = rdistphi - ...
															frdata.rphiadjval(freq_index);
			update_ui_str(handles.RSPLText, ...
									sprintf('%.4f', dbspl(mags{R}(freq_index, rep))));
			% store attenuation
			atten{R}(freq_index, rep) = Ratten;
			if cal.CheckCal == R
				[tmpmag, tmpphi] = fitsinvec(resp{inChan}(start_bin:end_bin), ...
																1, iodev.Fs, freq);
				mags{REF}(freq_index, rep) = VtoPa * RMSsin * tmpmag;
				phis{REF}(freq_index, rep) = tmpphi;
				[tmpdistmag, tmpdistphi] = ...
									fitsinvec(resp{inChan}(start_bin:end_bin), 1, ...
													iodev.Fs, 2*freq);
				dists{REF}(freq_index, rep) = tmpdistmag;
				distphis{REF}(freq_index, rep) = tmpdistphi;					
				fprintf('ref mag: %f    R mag: %f', ...
									dbspl(mags{REF}(freq_index, rep)), ...
									dbspl(mags{R}(freq_index, rep)) );
			elseif cal.CheckCal == BOTH
				[tmpmag, tmpphi] = fitsinvec(resp{inChan}(start_bin:end_bin), ...
															1, iodev.Fs, freq);
				mags{REFR}(freq_index, rep) = VtoPa * RMSsin * tmpmag;
				phis{REFR}(freq_index, rep) = tmpphi;
				[tmpdistmag, tmpphi] = ...
									fitsinvec(resp{inChan}(start_bin:end_bin), 1, ...
													iodev.Fs, 2*freq);
				dists{REFR}(freq_index, rep) = tmpdistmag;
				distphis{REFR}(freq_index, rep) = tmpdistphi;
				fprintf('refR mag: %f    R mag: %f', ...
									dbspl(mags{REFR}(freq_index, rep)), ...
									dbspl(mags{R}(freq_index, rep)) );
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
			% if MeasureLeak is requested, measure it!
			if handles.cal.MeasureLeak
				[lleakmag, lleakphi] = ...
									fitsinvec(resp{L}(start_bin:end_bin), 1, ...
														iodev.Fs, freq);
				[lleakdistmag, lleakdistphi] = ...
							fitsinvec(resp{L}(start_bin:end_bin), 1, ...
							iodev.Fs, 2*freq);
				% update L text display
				update_ui_str(handles.LValText, sprintf('%.4f', 1000*lleakmag));
				% compute distortion measures before applying corrections
				leakdists{L}(freq_index, rep) = lleakdistmag / lleakmag;
				% adjust for the gain of the preamp and apply correction
				% factors for RMS and microphone calibration
				lleakmag = RMSsin * lleakmag / ...
													(Gain(L)*frdata.lmagadjval(freq_index));
				% convert to Pascals (rms) and adjust phase measurements
				leakmags{L}(freq_index, rep) = VtoPa*(lleakmag);
				leakphis{L}(freq_index, rep) = lleakphi - ...
															frdata.lphiadjval(freq_index);
				leakdistphis{L}(freq_index, rep) = ...
											lleakdistphi - frdata.lphiadjval(freq_index);
				update_ui_str(handles.LSPLText, ...
										sprintf('%.4f', ...
												dbspl(leakmags{L}(freq_index, rep))));
			else
				update_ui_str(handles.LValText, '---');
				update_ui_str(handles.LSPLText, '---');				
			end
			% plot the response
			if handles.cal.MeasureLeak
				Lacq = downsample(resp{L}, handles.cal.deciFactor);
				refreshdata(H.Lacq, 'caller');
				[tmpf, Lfft] = daqdbfft(resp{L}(start_bin:end_bin), iodev.Fs, ...
												length(resp{L}(start_bin:end_bin))); %#ok<*ASGLU>
				refreshdata(H.Lfft, 'caller');
			end
			Racq = downsample(resp{R}, handles.cal.deciFactor);
			refreshdata(H.Racq, 'caller');

			[tmpf, Rfft] = daqdbfft(resp{R}(start_bin:end_bin), iodev.Fs, ...
												length(resp{R}(start_bin:end_bin)));
			refreshdata(H.Rfft, 'caller');
			drawnow
			% draw spectrogram
			axes(handles.Rspecgram);
			myspectrogram(resp{R}, iodev.Fs, ...
									[10 5], @hamming, handles.SpectrumWindow, ...
									[-100 -1], false, 'default', false, 'per');
			% save raw data
			if handles.cal.SaveRawData
				fp = fopen(rawfile, 'a');
				if handles.cal.MeasureLeak
					writeCell(fp, resp);
				else
					writeCell(fp, resp(R));
				end
				fclose(fp);
			end
			% pause for ISI (convert to seconds)
			pause(0.001*cal.ISI);
			% check for abort button press
			if read_ui_val(handles.AbortCtrl) == 1
				% if so, stop
				disp('abortion detected')
				break
			end
		end
		%---------------------------------------------------------------------
		% now, collect the background data for frequency FREQ, RIGHT channel
		%---------------------------------------------------------------------
		if handles.cal.CollectBackground
			Lstim = Nullstim_downsample;
			Rstim = Nullstim_downsample;
			refreshdata(H.Lstim, 'caller');
			refreshdata(H.Rstim, 'caller');
			for rep = 1:cal.Nreps
				% update the reps display value
				update_ui_str(handles.RepNumText, sprintf('%d R (bg)', rep));
				% play the sound;
				if handles.DAQSESSION
					[resp, indx] = handles.iofunction(iodev, Satt, ...
															handles.cal.SweepDuration);
				else
					[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
				end
				% filter the data if asked
				if handles.cal.InputFilter
					tmp = sin2array(resp{L}, 1, iodev.Fs);
					resp{L} = filtfilt(handles.cal.fcoeffb, ...
													handles.cal.fcoeffa, tmp);
					tmp = sin2array(resp{R}, 1, iodev.Fs);
					resp{R} = filtfilt(handles.cal.fcoeffb, ...
													handles.cal.fcoeffa, tmp);
					clear tmp
				end
				% determine the magnitude and phase of the response
				bgmag = fitsinvec(resp{inChan}(start_bin:end_bin), 1, ...
														iodev.Fs, freq);
				% update RValText
				update_ui_str(handles.RValText, sprintf('%.4f', 1000*bgmag));
				% adjust for the gain of the preamp and apply correction
				% factors for RMS and microphone calibration
				bgmag_adjusted = RMSsin * bgmag / ...
												(Gain(R)*frdata.rmagadjval(freq_index));
				% Store the values in the cell arrays for later averaging
				% (we'll do the averages later in order to save time while
				%  running the calibration curves)
				% adjust for the gain of the preamp and convert to Pascals
				bgmags{R}(freq_index, rep) = VtoPa*(bgmag_adjusted);
				update_ui_str(handles.RSPLText, ...
									sprintf('%.4f', dbspl(bgmags{R}(freq_index, rep))));
				update_ui_str(handles.LValText, '---');
				update_ui_str(handles.LSPLText, '---');
				% plot the response
				Lacq = downsample(resp{L}, handles.cal.deciFactor);
				Racq = downsample(resp{R}, handles.cal.deciFactor);
				refreshdata(H.Lacq, 'caller');
				refreshdata(H.Racq, 'caller');
				% plot fft
				figure(10)
				[tmpf, tmpm] = daqdbfft(resp{R}(start_bin:end_bin), ...
											iodev.Fs, length(resp{R}(start_bin:end_bin)));
				plot(tmpf, tmpm);
				title('Right Background')
				% save raw data if necessary
				if handles.cal.SaveRawData
					fp = fopen(rawfile, 'a');
					writeCell(fp, resp); 				
					fclose(fp);
				end
				% Pause for ISI
				pause(0.001*cal.ISI);
				% check for abort button press
				if read_ui_val(handles.AbortCtrl) == 1
					% if so, stop
					disp('abortion detected')
					break
				end
			end
		end
	% END OF R CAL
	end
	% check for abort button press
	if read_ui_val(handles.AbortCtrl) == 1
		% if so, stop
		disp('abortion detected')
		break
	end
	% increment frequency counter
	freq_index = freq_index + 1;
end %********************End of Cal loop

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Exit gracefully (close DAQ objects, etc)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_NIexit;

% % update the reps display value
% update_ui_str(handles.RepNumText, sprintf('%d R', rep));
%-------------------------------------------------------
% set COMPLETE if we made it to the last frequency, 
% otherwise, assume that ABORT was engaged and exit the run
%-------------------------------------------------------
if freq_index == Nfreqs+1
	COMPLETE = 1;
else
	COMPLETE = 0;
	return
end
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Compute Averages
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
freq_index = 1;
%*******************************LOOP through the frequencies
for F = 1:Nfreqs
	freq = Freqs(F);
	% compute the averages for this frequency

	% magnitude (dB) = db(rms) + atten
	magsraw{L}(freq_index, :) = dbspl(mags{L}(freq_index, :));
	magsraw{R}(freq_index, :) = dbspl(mags{R}(freq_index, :));
	mags{L}(freq_index, :) = dbspl(mags{L}(freq_index, :)) + ...
												atten{L}(freq_index, :);
	mags{R}(freq_index, :) = dbspl(mags{R}(freq_index, :)) + ...
												atten{R}(freq_index, :);
	% if Check data, save it
	if cal.CheckCal == L
		magsraw{REF}(freq_index, :) = dbspl(mags{REF}(freq_index, :));
		mags{REF}(freq_index, :) = dbspl(mags{REF}(freq_index, :)) + ...
														atten{L}(freq_index, :);
	elseif cal.CheckCal == R
		magsraw{REF}(freq_index, :) = dbspl(mags{REF}(freq_index, :));
		mags{REF}(freq_index, :) = dbspl(mags{REF}(freq_index, :)) + ...
														atten{R}(freq_index, :);
	elseif cal.CheckCal == BOTH
		magsraw{REFL}(freq_index, :) = dbspl(mags{REFL}(freq_index, :));
		magsraw{REFR}(freq_index, :) = dbspl(mags{REFR}(freq_index, :));
		mags{REFL}(freq_index, :) = dbspl(mags{REFL}(freq_index, :)) + ...
														atten{L}(freq_index, :);
		mags{REFR}(freq_index, :) = dbspl(mags{REFR}(freq_index, :)) + ...
														atten{R}(freq_index, :);
	end
	
	% store in caldata struct
	for channel = 1:handles.Nchannels				
		caldata.mag(channel, freq_index) = mean( mags{channel}(freq_index, :) );
		caldata.mag_stderr(channel, freq_index) = ...
													std( mags{channel}(freq_index, :) );

		caldata.phase(channel, freq_index) = ...
											mean( unwrap(phis{channel}(freq_index, :)) );
		caldata.phase_stderr(channel, freq_index) = ...
											std( unwrap(phis{channel}(freq_index, :)) );

		caldata.dist(channel, freq_index) = mean( dists{channel}(freq_index, :) );
		caldata.dist_stderr(channel, freq_index) = ...
														std( dists{channel}(freq_index, :) );
	end
	
	% store leak data if collected
	if handles.cal.MeasureLeak
		% compute the averages for this frequency
		leakmags{L}(freq_index, :) = dbspl(leakmags{L}(freq_index, :));
		leakmags{R}(freq_index, :) = dbspl(leakmags{R}(freq_index, :));
		leakphis{L}(freq_index, :) = leakphis{L}(freq_index, :) - ...
																phis{R}(freq_index, :);
		leakphis{R}(freq_index, :) = leakphis{R}(freq_index, :) - ...
																phis{L}(freq_index, :);
		for channel = 1:handles.Nchannels
			caldata.leakmag(channel, freq_index) = ...
											mean( leakmags{channel}(freq_index, :) );
			caldata.leakmag_stderr(channel, freq_index) = ...
											std( leakmags{channel}(freq_index, :) );
			caldata.leakphase(channel, freq_index) = ...
									mean( unwrap(leakphis{channel}(freq_index, :)) );
			caldata.leakphase_stderr(channel, freq_index) = ...
									std( unwrap(leakphis{channel}(freq_index, :)) );
			caldata.leakdist(channel, freq_index) = ...
									mean( leakdists{channel}(freq_index, :) );
			caldata.leakdist_stderr(channel, freq_index) = ...
									std( leakdists{channel}(freq_index, :) );
			caldata.leakdistphis(channel, freq_index) = ...
									mean( leakdistphis{channel}(freq_index, :) );
			caldata.leakdistphis_stderr(channel, freq_index) = ...
									std( leakdistphis{channel}(freq_index, :) );
		end
		caldata.leakmags = leakmags;
	end
	% if background data were collected, process and save them
	if handles.cal.CollectBackground
		bgdb{L}(freq_index, :) = dbspl(bgmags{L}(freq_index, :));
		bgdb{R}(freq_index, :) = dbspl(bgmags{R}(freq_index, :));		
		for channel = 1:handles.Nchannels
			caldata.background(channel, freq_index) =...
														mean( bgdb{channel}(freq_index, :) );
			caldata.background_stderr(channel, freq_index) =...
														std( bgdb{channel}(freq_index, :) );
		end
	end
	% go to next frequency!
	freq_index = freq_index + 1;
end
% assign to caldata struct
caldata.magsraw = magsraw;
caldata.magsV = magsV;
caldata.atten = atten;
caldata.calibration_settings = cal;
% store leak data if collected
if handles.cal.MeasureLeak
	caldata.leakmags = leakmags;
end
% save raw background data if collected
if handles.cal.CollectBackground
	caldata.bgmagsraw = bgmags;
end
if DEBUG
	caldata.magsdbug = magsdbug;
	caldata.phisdbug = phisdbug;
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% save handles and data
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
handles.caldata = caldata;
guidata(hObject, handles);
if cal.AutoSave
	disp(['Saving calibration data in ' handles.cal.calfile ' ...']);
	save(handles.cal.calfile, '-MAT', 'caldata');
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% plot the calibration data
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
try
	PlotCal(caldata);
catch errMsg
	save('NICal_RunCalibration.err', 'errMsg', '-MAT')
	errMsg.stack
end

disp('Finished.')

