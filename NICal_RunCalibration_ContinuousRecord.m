%--------------------------------------------------------------------------
% NICal_RunCalibration_ContinuousRecord.m
%--------------------------------------------------------------------------
% Collects Data for specified duration
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	7 Nov 2014 from NICal_RunCalibration,	SJS
%
% Revisions:
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
% need to create data file from calfile name
[fpath, fname, fext] = fileparts(handles.cal.calfile);
matfile = fullfile(fpath, [fname, '.mat']);
if exist(matfile, 'file')
	resp = uiyesno('title', 'Save File', 'string', ...
							'File exists! Overwrite?', 'default', 'No');
	if strcmpi(resp, 'No')
		[pathstr, fname, fext] = fileparts(matfile);
		[newname, newpath] = uiputfile('*.mat', ...
													'Save calibration data to file', ...
													fullfile(pathstr, [fname '_1' fext]));
		if isequal(newname, 0) || isequal(newpath, 0)
			return
		else
			matfile= fullfile(newpath, newname);
		end
	end
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
% write cal struct to rawdata file
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
if handles.cal.SaveRawData
	[pathstr, fname, fext] = fileparts(matfile);
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


%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% create null stimulus and time vector for plots, set up plots
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% fake acquired data
zeroacq = syn_null(handles.cal.SweepDuration, handles.iodev.Fs, 0);
zeroacq = downsample(zeroacq, handles.cal.deciFactor);
acqpts = length(zeroacq);
% time vector for stimulus plots
dt = 1/handles.iodev.Fs;
tvec_acq = 1000*dt*(0:(acqpts-1));
% compute # of points per sweep
SweepPoints = ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);

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
% acq
Lacq = zeroacq;
Racq = zeroacq;
% FFT
nfft = SweepPoints;
tmp = zeros(1, nfft);
[fvec, Lfft] = daqdbfft(tmp, handles.iodev.Fs, nfft);
[fvec, Rfft] = daqdbfft(tmp, handles.iodev.Fs, nfft);
% convert fvec to kHz
fvec = 0.001 * fvec;
clear tmp

%-------------------------------------------------------
% plot null data, save handles for time-domain plots
%-------------------------------------------------------
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

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Acquire data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
fprintf('Ready to acquire data...\n\n\n');

%----------------------------------------------------------------
% open a menu-panel for user to start (or cancel) acquisition
%----------------------------------------------------------------
userResp = menu('Acquiring 1 Sweep of Data', 'Start', 'Cancel');

%----------------------------------------------------------------
% if user hit cancel, stop
%----------------------------------------------------------------
if userResp == 2
	fprintf('Cancelling acquisition!\n\n');
	disp('...closing NI devices...');

	NICal_NIexit;
	COMPLETE = 0;
	return
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Now initiate sweep
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%---------------------------------------------------
% make a local copy of the cal settings structure
%---------------------------------------------------
cal = handles.cal;
%---------------------------------------------------
% make local copy of iodev  control struct
%---------------------------------------------------
iodev = handles.iodev;

		
%-------------------------------------------------------
% Build fake stimulus output array 
%-------------------------------------------------------
Satt = syn_null(cal.StimDuration, iodev.Fs, 2);

%-------------------------------------------------------
% now, collect the data 
%-------------------------------------------------------
inChan = 1;

% get data, (play null stim);
[resp, indx] = handles.iofunction(iodev, Satt, SweepPoints);

rawresp = resp;
% filter the data if asked
if handles.cal.InputFilter
	tmp = sin2array(resp{L}, 1, iodev.Fs);
	resp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
end
	
% plot the response and FFT
Lacq = downsample(resp{L}, handles.cal.deciFactor);
refreshdata(H.Lacq, 'caller');
[tmpf, Lfft] = daqdbfft(resp{L}, iodev.Fs, length(resp{L}));
refreshdata(H.Lfft, 'caller');
drawnow

% draw spectrogram
axes(handles.Lspecgram);
myspectrogram(resp{L}, iodev.Fs, ...
								[10 5], @hamming, handles.SpectrumWindow, ...
								[-100 -1], false, 'default', false, 'per');
			
% save raw data
if handles.cal.SaveRawData
	fp = fopen(rawfile, 'a');
	writeCell(fp, rawresp);
	fclose(fp);
end
			
% Pause for ISI
pause(0.001*cal.ISI);			
		
%---------------------------------------------------------------------
% now, collect the background data for frequency FREQ, LEFT channel
%---------------------------------------------------------------------
if handles.cal.CollectBackground
	% play the sound;
	[bgresp, indx] = handles.iofunction(iodev, Satt, SweepPoints);
	rawbgresp = bgresp;
	% filter the data if asked
	if handles.cal.InputFilter
		tmp = sin2array(bgresp{L}, 1, iodev.Fs);
		bgresp{L} = filtfilt(handles.cal.fcoeffb, handles.cal.fcoeffa, tmp);
		clear tmp
	end
	
	pause(0.001*cal.ISI);

	% save raw data
	if handles.cal.SaveRawData
		fp = fopen(rawfile, 'a');
		writeCell(fp, rawbgresp);
		fclose(fp);
	end
end


%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Exit gracefully (close TDT objects, etc)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
NICal_NIexit;

COMPLETE = 1;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% save the calibration data
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
if handles.cal.CollectBackground
	save(matfile, '-MAT', 'cal', 'resp', 'bgresp');
else
	save(matfile, '-MAT', 'cal', 'resp');	
end
disp('Finished.')

figure(1)
c = processContinuousData(matfile);
subplot(212)
ylim([40 1.1*max(c.dbvals)])
subplot(211)
title(matfile, 'Interpreter', 'none');


