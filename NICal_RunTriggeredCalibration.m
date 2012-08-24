%--------------------------------------------------------------------------
% NICal_RunTriggeredCalibration.m
%--------------------------------------------------------------------------
% Runs the TTL triggered speaker calibration
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	18 July, 2012 from NICal_RunCalibration (originally named
% 				NICal_TriggeredCalibration)
%
% Revisions:
%	24 Aug 2012 (SJS): renamed NICal_RunTriggeredCalibration to
% 							 more accurately reflect function as well as to
% 							 mirror NICal_RunCalibration
%--------------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Constants
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%---------------------------------------------
% Global
%---------------------------------------------
NICal_Constants;

%---------------------------------------------
% Local
%---------------------------------------------
% set the COMPLETE flag to 0
COMPLETE = 0;

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Initialization Scripts
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%---------------------------------------------
%---------------------------------------------
% Load Microphone calibration data
%---------------------------------------------
%---------------------------------------------
if read_ui_val(handles.FRenableCtrl) == 0
	DAscale = read_ui_str(handles.DAscaleCtrl, 'n');
	handles.cal.mic_fr = [];
	handles.cal.DAscale = DAscale;
else
	load(handles.cal.mic_fr_file, 'frdata');
	if ~isfield(frdata, 'DAscale')
		frdata.DAscale = frdata.calsettings.DAscale;
	end
	handles.cal.mic_fr = frdata;
end
% save the GUI handle information
guidata(hObject, handles);

%---------------------------------------------
%---------------------------------------------
% settings
%---------------------------------------------
%---------------------------------------------

% read in the gain on the mic preamp
Gain_dB = handles.cal.MicGain;
% convert dB to linear scale
Gain = invdb(Gain_dB);
% this is the sensitivity of the calibration mic in V / Pa
% if FR file is used, get sens. from there, 
if read_ui_val(handles.FRenableCtrl) == 1
	CalMic_sense = frdata.calsettings.CalMic_sense;
else
	% otherwise, use cal information
	CalMic_sense = handles.cal.MicSensitivity;
end
% pre-compute the V -> Pa conversion factor
VtoPa = (CalMic_sense^-1);
% precompute the volts -> RMS conversion factor for sinusoids (0.7071)
RMSsin = 1/sqrt(2);

%---------------------------------------------
% KLUDGE!!!!!!!
%---------------------------------------------
handles.Nchannels = 2;
guidata(hObject, handles);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% setup output files - automagic triggering will write data to
% binary output file
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
tmppath = fileparts(handles.cal.calfile);
tmpname = sprintf('NICaldata-%s-1.daq', date);
tmpfilename = fullfile(tmppath, tmpname);

[filename, pathname] = uiputfile(	{'*.daq', 'Matlab DAQ file (*.daq)'}, ...
												'Save Triggered data', ...
												tmpfilename );
											
if isequal(filename, 0) || isequal(pathname, 0)
	disp('Cancelling')
	return
else
	% parse filenames
	handles.OutputDataPath = pathname;
	handles.OutputDataFile = filename;
	[~, filebase] = fileparts(handles.OutputDataFile);
	filebase = filebase(1:(length(filebase)-2));
	handles.OutputMatFile = [filebase '.mat'];
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Start DAQ things
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
[handles, init_status] = NICal_NIinit_triggeredacq(handles);
guidata(hObject, handles);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Define a bandpass filter for processing the data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Nyquist frequency
fnyq = handles.cal.Fs/2;
% passband definition
handles.cal.fband = [handles.cal.InputHPFc handles.cal.InputLPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[handles.cal.fcoeffb, handles.cal.fcoeffa] = ...
					butter(handles.cal.forder, handles.cal.fband, 'bandpass');

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Setup caldata struct for storing the calibration data
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% NICal_caldata_init;


%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% create null response and time vector for plots, set up plots, global
% vars for callback
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
global tvec_acq Lacq Racq fvec Lfft Rfft H SweepPoints deciFactor start_bin end_bin

% plot decimation factor
deciFactor = handles.deciFactor;

% fake acquired data
zeroacq = syn_null(handles.cal.SweepDuration, handles.iodev.Fs, 0);
zeroacq = downsample(zeroacq, deciFactor);
acqpts = length(zeroacq);
% time vector for stimulus plots
tvec_acq = 1000*(1/handles.iodev.Fs)*(0:(acqpts-1));
% compute # of points per sweep
SweepPoints = ms2samples(handles.cal.SweepDuration, handles.iodev.Fs);

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Setup Plots
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%-------------------------------------------------------
% create arrays for plotting and plot them
%-------------------------------------------------------
% acq
Lacq = zeroacq;
Racq = zeroacq;
% FFT
nfft = length(start_bin:end_bin);
tmp = zeros(1, nfft);
[fvec, Lfft] = daqdbfft(tmp, handles.iodev.Fs, nfft);
[fvec, Rfft] = daqdbfft(tmp, handles.iodev.Fs, nfft);
% convert fvec to kHz
fvec = 0.001 * fvec;
clear tmp

%-------------------------------------------------------
% plot null data, save handles for time-domain plots
%-------------------------------------------------------
% stimulus
H.Lstim = plot(0);
H.Rstim = plot(0);
% response
H.Lacq = plot(handles.Lmicplot, tvec_acq, Lacq, 'g');
set(H.Lacq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Lacq');
xlabel(handles.Lmicplot, 'Time (ms)')
ylabel(handles.Lmicplot, 'V')
H.Racq = plot(handles.Rmicplot, tvec_acq, Racq, 'r');
set(H.Racq, 'XDataSource', 'tvec_acq', 'YDataSource', 'Racq');
xlabel(handles.Rmicplot, 'Time (ms)')

%-------------------------------------------------------
% plot null data, save handles for frequency-domain plots
%-------------------------------------------------------
H.Lfft = plot(handles.Lfftplot, fvec, Lfft);
set(H.Lfft, 'XDataSource', 'fvec', 'YDataSource', 'Lfft');
xlabel(handles.Lfftplot, 'Frequency (kHz)')
ylabel(handles.Lfftplot, 'dBV')
H.Rfft = plot(handles.Rfftplot, fvec, Rfft);
set(H.Rfft, 'XDataSource', 'fvec', 'YDataSource', 'Rfft');
xlabel(handles.Rfftplot, 'Frequency (kHz)')

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% set the start and end bins for the calibration
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
start_bin = ms2bin(handles.cal.StimDelay + handles.cal.StimRamp, handles.iodev.Fs);
if ~start_bin
	start_bin = 1;
end
end_bin = start_bin + ms2bin(handles.cal.StimDuration-handles.cal.StimRamp, handles.iodev.Fs);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Acquire data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
fprintf('Ready to acquire data...\n\n\n');

stopFlag = 0;
rep = 1;
freq_index = 1;

userResp = menu('Acquiring Triggered Data', 'Start', 'Cancel');

if userResp ~= 2
	runFLAG = 1;
	while runFLAG
		%START ACQUIRING
		start(handles.iodev.NI.ai);
		fprintf('Writing to file %s \n', handles.iodev.NI.ai.LogFileName);

		try
			wait(handles.iodev.NI.ai, handles.cal.TriggerTimeout);
		catch errEvent
			disp('TIMEOUT!');
			runFLAG = 0;
		end
	end
	% stop acquiring
	fprintf('... terminating\n')
	stop(handles.iodev.NI.ai);
end

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Exit gracefully (close TDT objects, etc)
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Close hardware
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
disp('...closing NI devices...');

% get event log
EventLogAI = showdaqevents(handles.iodev.NI.ai);

% delete and clear ai and ch0 object
delete(handles.iodev.NI.ai);
delete(handles.iodev.NI.chI);
clear handles.iodev.NI.ai handles.iodev.NI.chI

% save settings information to mat file
save(fullfile(pwd, 'NICal_EventLogs.mat'), ...
		'EventLogAI'			, ...
		'-MAT' );

% save settings information to mat file
Fs = handles.iodev.Fs;
InputChannel = handles.cal.InputChannel;
MicChannel = handles.cal.Side;
save(fullfile(handles.OutputDataPath, handles.OutputMatFile), ...
		'InputChannel'			, ...
		'MicChannel'			, ...
		'acqpts'					, ...
		'VtoPa'					, ...
		'Fs'						, ...
		'EventLogAI'			, ...
		'-MAT' );

COMPLETE = 1;

disp('Finished.')


