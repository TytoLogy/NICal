%--------------------------------------------------------------------------
% NICal_MicrophoneCalibrate.m
%--------------------------------------------------------------------------
%  Script that runs the microphone calibration protocol
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	16 August, 2012 (SJS) from MicrophoneCal_RunSequential
%
% Revisions:
%
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
% get output  file - if it exists, check with user
%-----------------------------------------------------------------------
[pathstr, fname, fext] = fileparts(handles.cal.calfile);
% store original file to restore to uictrl after frdata run
calfile_orig = handles.cal.calfile;

[newname, newpath] = uiputfile('*.fr','Save microphone freq response data to file', pathstr);
if isequal(newname, 0) || isequal(newpath, 0)
	fprintf('Cancelling .fr file collection');
	return
else
	handles.cal.calfile = fullfile(newpath, newname);
	update_ui_str(handles.CalFileCtrl, handles.cal.calfile);
	guidata(hObject, handles);
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Initialization Scripts
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Load the settings and constants for NICal 
%--------------------------------------------------------------------------
NICal_settings;
%---------------------------------------------
% KLUDGE!!!!!!!
%---------------------------------------------
handles.cal.Nchannels = 2;
guidata(hObject, handles);

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
try
	[handles.cal.fcoeffb, handles.cal.fcoeffa] = ...
						butter(handles.cal.forder, handles.cal.fband, 'bandpass');
catch errMsg
	fprintf('handles.cal.fband = %f\n', handles.cal.fband);
	NICal_NIexit;
	errMsg
	return
end

guidata(hObject, handles);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% setup storage variables
%------------------------------------------------------------------------
%------------------------------------------------------------------------
NICal_frdata_init;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Read in the BK mic xfer function for pressure field 
% and get correction values for use with the free field mic
% If free-field, set correction factor to 1
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if handles.cal.FieldType == 2
	% interpolate to get values at desired freqs (data in dB)
	frdata.bkpressureadj = interp1(	handles.cal.bkdata.Response(:, 1), ...
												handles.cal.bkdata.Response(:, 2), ...
												Freqs);
	% convert to factor
	frdata.bkpressureadj = 1./invdb(frdata.bkpressureadj);
else
	frdata.bkpressureadj = ones(size(Freqs));
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set the start and end bins for the calibration
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
start_bin = ms2bin(handles.cal.StimDelay + handles.cal.StimRamp, handles.iodev.Fs);
if ~start_bin
	start_bin = 1;
end
end_bin = start_bin + ms2bin(handles.cal.StimDuration - handles.cal.StimRamp, handles.iodev.Fs);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% create null stimulus and time vector for plots, set up plots
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% create null stimulus
zerostim = syn_null(handles.cal.StimDuration, handles.iodev.Fs, 0);
% insert stim delay
zerostim = insert_delay(zerostim, handles.cal.StimDelay, handles.iodev.Fs);
% downsample (no need to plot all points)
zerostim = downsample(zerostim, deciFactor);
% downsample-factor adjusted sample interval
dt = deciFactor/handles.iodev.Fs;
% # output points
outpts = length(zerostim);
% time vector for stimulus plots
tvec_stim = 1000*dt*(0:(outpts-1));
% stimulus start and end points
stim_start = ms2bin(handles.cal.StimDelay, handles.iodev.Fs);
stim_end = stim_start + outpts - 1;
% fake acquired data
zeroacq = syn_null(handles.cal.SweepDuration, handles.iodev.Fs, 0);
zeroacq = downsample(zeroacq, deciFactor);
acqpts = length(zeroacq);
% time vector for stimulus plots
tvec_acq = 1000*dt*(0:(acqpts-1));
% compute # of points per sweep
SweepPoints = ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Build null output array
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% synthesize the L sine wave;
Nullstim = syn_null(handles.cal.StimDuration, handles.iodev.Fs, 1);
% scale the sound
Nullstim = 0 * Nullstim;
% insert delay
Nullstim = insert_delay(Nullstim, handles.cal.StimDelay, handles.iodev.Fs);
Nullstim_downsample =  downsample(Nullstim(1, :), deciFactor);

%-------------------------------------------------------
% create arrays for plotting and plot them
%-------------------------------------------------------
% stim
Lstim = zerostim;
Rstim = zerostim;
% acq
Lacq = zeroacq;
Racq = zeroacq;
H.Lstim = plot(handles.Lstimplot, tvec_stim, Lstim, 'g');
set(H.Lstim, 'XDataSource', 'tvec_stim', 'YDataSource', 'Lstim');
H.Rstim = plot(handles.Rstimplot, tvec_stim, Rstim, 'r');
set(H.Rstim, 'XDataSource', 'tvec_stim', 'YDataSource', 'Rstim');
H.Lacq = plot(handles.Lmicplot, tvec_acq, Lacq, 'g');
set(H.Lacq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Lacq');
H.Racq = plot(handles.Rmicplot, tvec_acq, Racq, 'r');
set(H.Racq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Racq');
set(handles.Lstimplot, 'XTickLabel', '');
set(handles.Lmicplot, 'XTickLabel', '');

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% pre-allocate some data cell arrays
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
mags = cell(1, handles.cal.Nchannels);
for c = 1:handles.cal.Nchannels
	mags{c} = zeros(handles.cal.Nfreqs, handles.cal.Nreps);
end
phis = mags; 
dists = mags;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% setup attenuation
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
if handles.cal.AttenFix 
	Latten = handles.cal.AttenFixValue;
	Ratten = handles.cal.AttenFixValue;
else
	warning('NICal:Atten', [mfilename ': AttenFix not set, using default StartAtten value']);
	% set the adjustable starting attenuator values
	Latten = handles.cal.StartAtten;
	Ratten = handles.cal.StartAtten;
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% another kludge... use channel 2 (RIGHT) as reference
%-----------------------------------------------------------------------
if handles.cal.Nchannels == 2
	REFCHAN = 2;
elseif handles.cal.Nchannels == 3
	REFCHAN = 3;
else
	REFCHAN = 2;
end
%-----------------------------------------------------------------------

%---------------------------------------------------
% make a local copy of the cal settings structure
%---------------------------------------------------
cal = handles.cal;
%---------------------------------------------------
% make local copy of iodev  control struct
%---------------------------------------------------
iodev = handles.iodev;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Now initiate sweeps
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% First, get the background noise level
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% plot the stim arrays
Lstim = zerostim;
Rstim = zerostim;
refreshdata(H.Lstim, 'caller');
refreshdata(H.Rstim, 'caller');

% pre-allocate the background data cell array
background = cell(1, handles.cal.Nchannels);
for c = 1:handles.cal.Nchannels
	background{c} = zeros(1, handles.cal.Nreps);
end
	
% pause to let things settle down
disp([mfilename ': collecting background data']);
update_ui_str(handles.FreqValText, 'Backgnd');
pause(1);

% check for abort button press
if read_ui_val(handles.AbortCtrl) == 1
	% if so, stop
	disp('abortion detected')
	break
end

for rep = 1:handles.cal.Nreps
	% update the reps display value
	update_ui_str(handles.RepNumText, sprintf('%d L', rep));

	% play the "sound"
	[resp, indx] = handles.iofunction(iodev, Nullstim, SweepPoints);

	% filter the data if asked
	if handles.cal.InputFilter
		tmp = sin2array(resp{L}, FILTSMOOTH_MS, iodev.Fs);
		resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
		tmp = sin2array(resp{R}, FILTSMOOTH_MS, iodev.Fs);
		resp{R} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
		clear tmp
	end
		
	% plot the response
	Lacq = downsample(resp{L}, deciFactor);
	Racq = downsample(resp{R}, deciFactor);
	refreshdata(H.Lacq, 'caller');
	refreshdata(H.Racq, 'caller');

	% determine the magnitude of the response
	% first, the test channel (L)
	background{L}(rep) = rms(resp{L}(start_bin:end_bin) / Gain(L));
	background{REFCHAN}(rep) = rms(resp{REFCHAN}(start_bin:end_bin) / Gain(REFCHAN));
	% compute dB SPL using reference channel (L)
	background{REFCHAN}(rep) = dbspl(VtoPa*background{REFCHAN}(rep));

	% update display text, rms values
	update_ui_str(handles.LValText, sprintf('%.4f', 1000*background{L}(rep)));
	update_ui_str(handles.RValText, sprintf('%.4f', 1000*background{REFCHAN}(rep)));
	% update display text, SPL values
	update_ui_str(handles.LSPLText, '---');
	update_ui_str(handles.RSPLText, sprintf('%.2f', background{REFCHAN}(rep)));

	% store the response data
	rawdata.background{rep} = cell2mat(resp');

	% Pause for ISI
	pause(0.001*cal.ISI);
end

% average and st. dev
for c = 1:handles.cal.Nchannels
	frdata.background(c, 1) = mean( background{c} );
	frdata.background(c, 2) = std( background{c} );
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%*******************************LOOP through the frequencies
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
	
% PAUSE	
disp([mfilename ': now running calibration...']);
pause(1)
	
tic
handles.STOP_FLG = 0;	
rep = 1;
freq_index = 1;

while ~handles.STOP_FLG && freq_index <= handles.cal.Nfreqs
	% get the current frequency
	freq = Freqs(freq_index);

	% update the frequency display value
	update_ui_str(handles.FreqValText, sprintf('%d', round(freq)));

	%-------------------------------------------------------
	% Build stimulus output array for this frequency
	%-------------------------------------------------------
	% synthesize and scale the Left output channel sine wave;
	[S, stimspec.RMS, stimspec.phi] = syn_calibrationtone(cal.StimDuration, iodev.Fs, freq, 0, 'L');
	S = cal.DAscale*S;
	% apply the sin^2 amplitude envelope eliminate
	% transients in the stimulus
	S = sin2array(S, cal.StimRamp, iodev.Fs);
	% insert delay
	S = insert_delay(S, cal.StimDelay, iodev.Fs);
	% save in Satt
	Satt = S;
	% attenuate signale
	Satt(L, :) = handles.attfunction(S(L, :), Latten);
	Satt(R, :) = handles.attfunction(S(R, :), MAX_ATTEN);
	update_ui_str(handles.LAttenText, Latten);
	update_ui_str(handles.RAttenText, MAX_ATTEN);
	% plot the stimuli - set R stim to zero
	Lstim = downsample(Satt(L, :), deciFactor);
	Rstim = downsample(Satt(R, :), deciFactor);
	refreshdata(H.Lstim, 'caller');
	refreshdata(H.Rstim, 'caller');
	
	%-------------------------------------------------------
	% now, collect the data for frequency FREQ, LEFT channel
	%-------------------------------------------------------
	for rep = 1:handles.cal.Nreps
		% update the reps display value
		update_ui_str(handles.RepNumText, sprintf('%d', rep));

		% play the sound;
		[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
		
		% filter the data if asked
		if handles.cal.InputFilter
			tmp = sin2array(resp{L}, FILTSMOOTH_MS, iodev.Fs);
			resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
			tmp = sin2array(resp{R}, FILTSMOOTH_MS, iodev.Fs);
			resp{R} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
			clear tmp
		end

		% determine the magnitude and phase of the response
		[lmag, lphi] = fitsinvec(resp{L}(start_bin:end_bin), 1, iodev.Fs, freq);
		[rmag, rphi] = fitsinvec(resp{R}(start_bin:end_bin), 1, iodev.Fs, freq);
		[refmag, refphi] = fitsinvec(resp{REFCHAN}(start_bin:end_bin), 1, iodev.Fs, freq);
		[ldistmag, ldistphi] = fitsinvec(resp{L}(start_bin:end_bin), 1, iodev.Fs, 2*freq);				
		[rdistmag, rdistphi] = fitsinvec(resp{R}(start_bin:end_bin), 1, iodev.Fs, 2*freq);				
		[refdistmag, refdistphi] = fitsinvec(resp{REFCHAN}(start_bin:end_bin), 1, iodev.Fs, 2*freq);			
		update_ui_str(handles.LValText, sprintf('%.4f', 1000*lmag));
		update_ui_str(handles.RValText, sprintf('%.4f', 1000*refmag));

		% compute 2nd harmonic distortion ratio
		dists{L}(freq_index, rep) = ldistmag / lmag;
		dists{R}(freq_index, rep) = rdistmag / rmag;
		dists{REFCHAN}(freq_index, rep) = refdistmag / refmag;

		% adjust for the gain of the preamp and convert to RMS
		lmag = RMSsin * lmag / Gain(L);
		rmag = RMSsin * rmag / Gain(R);
		% adjust REFerence mic for gain, adjust for free-field or 
		% pressure-field, and convert to RMS
		refmag = RMSsin * refmag * frdata.bkpressureadj(freq_index) ./ Gain(REFCHAN);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% DEBUGGING
		% frdata.refmag(freq_index, rep) = refmag;
		% frdata.RMSsin(freq_index, rep) = RMSsin;
		% frdata.GainRef(freq_index, rep) = Gain(REFCHAN);
		% frdata.adjfactor(freq_index, rep) =  frdata.bkpressureadj(freq_index);
		% frdata.rawmags(freq_index, rep) =  RMSsin * refmag / Gain(REFCHAN);
		% frdata.adjmags(freq_index, rep) = RMSsin * frdata.bkpressureadj(freq_index) * refmag  ./ Gain(REFCHAN);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% store the data in arrays
		mags{L}(freq_index, rep) = lmag;
		mags{R}(freq_index, rep) = rmag;
		mags{REFCHAN}(freq_index, rep) = refmag;
		phis{L}(freq_index, rep) = lphi;
		phis{R}(freq_index, rep) = rphi;
		phis{REFCHAN}(freq_index, rep) = refphi;

		% plot the response
		Lacq = downsample(resp{L}, deciFactor);
		Racq = downsample(resp{REFCHAN}, deciFactor);
		refreshdata(H.Lacq, 'caller');
		refreshdata(H.Racq, 'caller');
	
		update_ui_str(handles.RSPLText, sprintf('%.4f', dbspl(VtoPa*mags{REFCHAN}(freq_index, rep))));
		
		% Check for possible clipping (values > 5V for NI)
		for channel = 1:handles.cal.Nchannels
			if max(abs(resp{channel})) >= CLIPVAL
				handles.STOP_FLG = channel;
			end
		end

		% store the raw response data
		rawdata.resp{freq_index, rep} = cell2mat(resp');

		% Pause for ISI
		pause(0.001*cal.ISI);
	end	%******************** END OF REPS

	% compute the averages for this frequency
	for channel = 1:handles.cal.Nchannels
		frdata.mag(channel, freq_index) = mean( mags{channel}(freq_index, :) );
		frdata.mag_stderr(channel, freq_index) = std( mags{channel}(freq_index, :) );
		frdata.phase(channel, freq_index) = mean( unwrap(phis{channel}(freq_index, :)) );
		frdata.phase_stderr(channel, freq_index) = std( unwrap(phis{channel}(freq_index, :)) );
		frdata.dist(channel, freq_index) = mean( dists{channel}(freq_index, :) );
	end

	% increment frequency index counter
	freq_index = freq_index + 1;

	% check if user pressed ABORT button or if STOP_FLG is set
	if handles.STOP_FLG
		disp('STOP_FLG detected')
		cal.timer = toc;
		break;
	end
	if read_ui_val(handles.AbortCtrl) == 1
		disp('abortion detected')
		cal.timer = toc;
		break
	end
end %********************End of Cal loop

% get the time
cal.timer = toc;

if handles.STOP_FLG
	errstr = sprintf('Possible clip on channel %d', handles.STOP_FLG);
	errordlg(errstr, 'Clip alert!');
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Exit gracefully (close TDT objects, etc)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_NIexit;

%-------------------------------------------------------
% set COMPLETE if we made it to the last frequency, 
% otherwise, assume that ABORT was engaged and exit the run
%-------------------------------------------------------
if freq_index == handles.cal.Nfreqs+1
	COMPLETE = 1;
else
	% if not, incomplete, skip the calculations and
	% return
	COMPLETE = 0;
	handles.frdata = frdata;
	handles.rawdata = rawdata;
	return
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% data are complete, do some computations
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% microphone adjust factors are:
% 	magnitude adj = knowles mic Vrms / Ref mic Vrms
% 	phase adj = Knowls mic deg - ref mic radians

% compute L and R magnitude correction
frdata.ladjmag = frdata.mag(R, :) ./ frdata.mag(REFCHAN, :);
frdata.radjmag = frdata.mag(R, :) ./ frdata.mag(REFCHAN, :);

% compute L and R phase correction
frdata.ladjphi = frdata.phase(R, :) - frdata.phase(REFCHAN, :);
frdata.radjphi = frdata.phase(R, :) - frdata.phase(REFCHAN, :);
	
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% save handles and data and temp file
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
handles.frdata = frdata;
handles.rawdata = rawdata;
handles.cal = cal;
guidata(hObject, handles);

save(handles.cal.calfile, 'frdata', 'cal', '-mat');
	
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% plot curves
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% average data (non-normalized)
figure
subplot(211)
hold on
	errorbar(frdata.freq, 1000*frdata.mag(L, :), 1000*frdata.mag_stderr(L, :), 'g');
	errorbar(frdata.freq, frdata.mag(REFCHAN, :), frdata.mag_stderr(REFCHAN, :), 'k-.');
hold off
ylabel('Magnitude')
legend('Left (X1000)', 'Ref');
title('Calibration Results')
subplot(212)
hold on
	errorbar(frdata.freq, frdata.phase(L, :), frdata.phase_stderr(L, :), 'g');
	errorbar(frdata.freq, frdata.phase(REFCHAN, :), frdata.phase_stderr(REFCHAN, :), 'k-.');
hold off
ylabel('Phase')
xlabel('Frequency (Hz)')
legend('Left', 'Ref');

% Normalized data plot	
figure
subplot(211)
plot(frdata.freq, normalize(frdata.mag(1, :)), 'g.-')
hold on
	plot(frdata.freq, normalize(frdata.mag(REFCHAN, :)), 'k.-')
hold off
ylabel('Normalized Magnitude')
legend('Left', 'Ref');
title('Normalized Frequency Response')
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'Color', .5*[1 1 1]);

subplot(212)
plot(frdata.freq, unwrap(frdata.phase(1, :)), 'g.-');
hold on
	plot(frdata.freq, unwrap(frdata.phase(REFCHAN, :)), 'k.-');
hold off
ylabel('Unwrapped Phase')
xlabel('Frequency (Hz)')
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(gca, 'Color', .5*[1 1 1]);
	
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% plot the calibration data
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
try
	PlotFR(frdata)
catch errMsg
	save('NICal_RunCalibration.err', 'errMsg', '-MAT')
	errMsg.stack
end

disp('Finished.')

