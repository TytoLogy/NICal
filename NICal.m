function varargout = NICal(varargin)
% NICAL M-file for NICal.fig
%      NICAL, by itself, creates a new NICAL or raises the existing
%      singleton*.
%
%      H = NICAL returns the handle to a new NICAL or the handle to
%      the existing singleton*.
%
%      NICAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NICAL.M with the given input arguments.
%
%      NICAL('Property','Value',...) creates a new NICAL or raises the
%      existing singleton*.  
%

%-------------------------------------------------------------------------
% Calibration Algorithm:
% 
% First, the sound level from the microphone is checked.  If the level 
% (measured in dB SPL) is too low, the value of the attenuator is
% decreased.  If it is too high, attenuation is increased.
% 
% Then, the calibration tone at given frequency is played from the
% speakers.  The tone is scaled to the maximum output level of the D/A 
% converter.
% 
% A phased lock loop (fitsinvec() function) is used to determine the 
% magnitude and phase of the response measured by the microphone.
% 
% The magnitude is then divided by the scaling factor
% from the microphone calibration.  This adjusts for the discrepancy
% between the earphone microphone (Knowles) and the reference
% (Bruel & Kjaer) microphone which is the reference, flat response.
% Phase is adjusted by subtracting the reference phase from the measured
% phase.
% 
% The magnitude is then converted to the RMS value 
% 
%-------------------------------------------------------------------------

% Last Modified by GUIDE v2.5 06-Feb-2017 15:51:45

% Begin initialization code - DO NOT EDIT
	gui_Singleton = 1;
	gui_State = struct('gui_Name',       mfilename, ...
					   'gui_Singleton',  gui_Singleton, ...
					   'gui_OpeningFcn', @NICal_OpeningFcn, ...
					   'gui_OutputFcn',  @NICal_OutputFcn, ...
					   'gui_LayoutFcn',  [] , ...
					   'gui_Callback',   []);
	if nargin && ischar(varargin{1})
		gui_State.gui_Callback = str2func(varargin{1});
	end

	if nargout
		[varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
	else
		gui_mainfcn(gui_State, varargin{:});
	end
% End initialization code - DO NOT EDIT
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% --- Executes just before NICal is made visible.
% Performs Initial Setup
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function NICal_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
	%----------------------------------------------------------
	%----------------------------------------------------------
	% Setup Paths
	%----------------------------------------------------------
	%----------------------------------------------------------
	disp([mfilename ': checking paths'])
	% behavior varies with platform
	if ispc
		% directory when using installed version:
		pdir = ['C:\TytoLogy\Toolboxes\TytoLogySettings\' getenv('USERNAME')];
	elseif ismac
		pdir = ['~/Work/Code/Matlab/dev/TytoLogy/TytoLogySettings/' ...
							getenv('USER')];
	end
	if isempty(which('ms2samples'))
		run(fullfile(pdir, 'tytopaths'))
	else
		disp([mfilename ': paths ok, launching programn'])
	end
	%----------------------------------------------------------
	%----------------------------------------------------------
	% check for version compatability
	%----------------------------------------------------------
	%----------------------------------------------------------
	tmp = textscan(version, '%s');
	tmp = tmp{1};
	handles.MATversion = tmp{1};
	if str2double(handles.MATversion(1)) < 8
		warning('Using legacy interface');
		handles.DAQSESSION = 0;
	else
		fprintf('Using Session interface');
		handles.DAQSESSION = 1;
	end
	guidata(hObject, handles);
	%----------------------------------------------------------
	%----------------------------------------------------------
	% load the configuration information and store as config structure
	% The HPSearch_Configuration.m function file will usually live in:
	%	<tytology base path>\<Toolboxes>TytoSettings\<username>\ 
	% directory
	%----------------------------------------------------------
	%----------------------------------------------------------
	% define user configuration path
	userconfigpath = ['C:\TytoLogy\Toolboxes\TytoLogySettings\' ...
								getenv('USERNAME') '\NICal\'];
	% load the configuration information, store in config structure
	if isempty(which('NICal_Configuration'))
		if ~exist(userconfigpath, 'dir')
			% if config isn't found, see if NICal should create it
			qstr = sprintf('User config directory %s does not exist!', ...
										userconfigpath);
			uresp = questdlg(	{qstr, 'Shall I create it?'}, ...
									'Configuration Not Found', ...
									'Yes', 'No', ...
									'Yes' );
			switch uresp
				case 'Yes'
					% create directory and create config from the program
					% default configuration
					mkdir(userconfigpath);
					addpath(userconfigpath);
					disp('Loading Defaults...');
					config = NICal_DefaultConfiguration;
				case 'No'
					% just load default
					disp('Loading Defaults...');
					config = NICal_DefaultConfiguration;
			end
		else
			addpath(userconfigpath);
			config = NICal_Configuration;
		end
	else
		config = NICal_Configuration;
	end
	% store configuration in handles;
	handles.config = config;
	handles.userconfigpath = userconfigpath;
	% save handles
	guidata(hObject, handles);	
	%----------------------------------------------------------
	%----------------------------------------------------------
	% Initial Calibration settings
	%----------------------------------------------------------
	%----------------------------------------------------------
	%---------------------------------------------
	% KLUDGE!!!!!!!
	%---------------------------------------------
	handles.Nchannels = 2;
	guidata(hObject, handles);
	%------------------------------------------------------------
	% first check to see if defaults file exists
	%------------------------------------------------------------
	INITDEFAULT = 0;
	if ~isempty(varargin)
		if any(strcmpi('InitDefaults', varargin))
			warndlg('Initializing Defaults!', mfilename);
			INITDEFAULT = 1;
		end
	end
	% load defaults
	handles.defaultsfile = fullfile(userconfigpath, [mfilename '_Defaults.mat']);
	if exist(handles.defaultsfile, 'file') && (INITDEFAULT == 0)
		cal = [];
		fprintf('Loading cal settings from defaults file %s ...\n', ...
							handles.defaultsfile)
		load(handles.defaultsfile);
		handles.cal = cal;
	else
		cal = NICal_calstruct_init;
		if INITDEFAULT
			save(handles.defaultsfile, 'cal');
		else
			% no defaults found, so use internal values
			disp('no defaults found, using internal values')
		end
		handles.cal = cal;
		clear cal;
	end
	guidata(hObject, handles);
	%----------------------------------------------------------
	% update user interface
	%----------------------------------------------------------
	NICal_UpdateUIFromCal(handles, handles.cal);
	%----------------------------------------------------------
	%----------------------------------------------------------
	% setup iodev NI device interface struct
	% using information from the config struct 
	%----------------------------------------------------------
	%----------------------------------------------------------
	% Dnum = device number - this is for dev1
	iodev.Dnum = sprintf('Dev%d', config.IODEVNUM);
	% reference channel
	iodev.REF = 1;
	iodev.status = 0;
	handles.iodev = iodev;
	guidata(hObject, handles);
	%----------------------------------------------------------
	%----------------------------------------------------------
	% set function handles from configuration data
	%----------------------------------------------------------
	%----------------------------------------------------------
	if handles.cal.TriggeredAcquisition == 1
		handles.initfunction = @nidaq_triggeredacq_init;
	else
		handles.initfunction = @nidaq_aiao_init;
	end
	handles.iofunction = config.IOFUNCTION;
	handles.attfunction = config.ATTENFUNCTION;
	guidata(hObject, handles);
	%----------------------------------------------------------
	% set default output path/file for calibration data
	%----------------------------------------------------------
	handles.cal.calfile = fullfile(config.DEFAULT_OUTPUT_PATH, 'nicaldata.cal');
	update_ui_str(handles.CalFileCtrl, handles.cal.calfile);
	guidata(hObject, handles);
	%--------------------------------------------------
	% spectrum settings
	%--------------------------------------------------
	handles.SpectrumWindow = 512;
	handles.ColorMap = 'gray';
	guidata(hObject, handles);
	%--------------------------------------------------
	% ToneStack settings
	%--------------------------------------------------
	handles.ToneStack = 0;
	set(handles.Menu_ToneStack, 'Checked', 'off');
	guidata(hObject, handles);
	%--------------------------------------------------
	% ToneSweep settings
	%--------------------------------------------------	
	set(handles.Menu_ToneSweep, 'Checked', 'off');
	handles.ToneSweep = 0;
	guidata(hObject, handles);
	%--------------------------------------------------
	% 	ContinuousRecord settings
	%--------------------------------------------------	
	set(handles.Menu_ContinuousRecord, 'Checked', 'off');
	handles.ContinuousRecord = 0;
	guidata(hObject, handles);
	%----------------------------------------------------------
	% Update handles structure
	%----------------------------------------------------------
	handles.CalComplete = 0;
	handles.output = hObject;
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% BUTTON CONTROL CALLBACKS
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Run Calibration callback
%-------------------------------------------------------------------------
function RunCalibrationCtrl_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
	%---------------------------------------------------------------
	% turn off calibration ctrl, enable abort ctrl
	%---------------------------------------------------------------
	disable_ui(handles.RunCalibrationCtrl);
	show_uictrl(handles.AbortCtrl);
	set(handles.AbortCtrl, 'Value', 0);
	%---------------------------------------------------------------
	% initialize complete flag - COMPLETE is used in the scripts
	% that do the heavy lifting...
	%---------------------------------------------------------------
	handles.CalComplete = 0;
	COMPLETE = 0; %#ok<NASGU>
	guidata(hObject, handles);
	%---------------------------------------------------------------
	% run appropriate calibration script
 	%---------------------------------------------------------------
	switch read_ui_val(handles.TriggeredAcquisitionCtrl) 
		case 0
			if handles.ToneStack
				NICal_RunCalibration_ToneStack
			elseif handles.ToneSweep
				NICal_RunCalibration_ToneSweep
			elseif handles.ContinuousRecord
				NICal_RunCalibration_ContinuousRecord
			else
				NICal_RunCalibration
			end
		case 1
			if handles.DAQSESSION
				NICal_RunTriggeredCalibration_Session;
			else
				NICal_RunTriggeredCalibration;
			end
	end
	%---------------------------------------------------------------
	% enable Calibration ctrl, disable abort ctrl
	%---------------------------------------------------------------
	enable_ui(handles.RunCalibrationCtrl);
	hide_uictrl(handles.AbortCtrl);
	set(handles.AbortCtrl, 'Value', 0);
	%---------------------------------------------------------------
	% save handles
	%---------------------------------------------------------------
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Abort running calibration
%--------------------------------------------------------------------------
function AbortCtrl_Callback(hObject, eventdata, handles) 
	disp('ABORTING Calibration!')
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% --- Executes on button press in MonitorCtrl.
%--------------------------------------------------------------------------
function MonitorCtrl_Callback(hObject, eventdata, handles)
	NICal_Monitor;
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Calibration Settings callbacks
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
% --- Executes on selection change in SideCtrl.
%--------------------------------------------------------------------------
function SideCtrl_Callback(hObject, eventdata, handles)
	handles.cal.Side = read_ui_val(hObject);
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------
% --- Executes on button press in TriggeredAcquisitionCtrl.
%-------------------------------------------------------------------------
function TriggeredAcquisitionCtrl_Callback(hObject, eventdata, handles)
	handles.cal.TriggeredAcquisition = read_ui_val(hObject);
	guidata(hObject, handles);
	% load list of handles to disable
	trig_handles_list;
	% need to select proper nidaq init function
	if handles.cal.TriggeredAcquisition == 1
		handles.initfunction = @nidaq_ai_init;
		disable_ui(trigger_handles);
	else
		handles.initfunction = @nidaq_aiao_init;
		enable_ui(trigger_handles);
	end
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% --- Executes on button press in FreqListCtrl.
%-------------------------------------------------------------------------
function FreqListCtrl_Callback(hObject, eventdata, handles)
	%------------------------------------------------------
	% get value of ui control (i.e., if checked or not)
	%------------------------------------------------------
	newVal = read_ui_val(hObject);
	%------------------------------------------------------
	% act accordingly
	%------------------------------------------------------
	if newVal
		% get filename and path to frequency list file (.m or .txt)
		[filename, pathname] = uigetfile(	...
														{	...
															'*.txt', 'text files (*.txt)'; ...
															'*.m', 'Matlab m-files (*.m)'; ...
															'*.*', 'All Files (*.*)'	...
														}, ...
														'Pick a file', ...
														handles.userconfigpath	...
													);
		% return if user pressed "Cancel" button
		if any([(filename == 0), (pathname == 0)])
			% uncheck box, set UseFreqList in handles to 0
			update_ui_val(hObject, 0);
			handles.cal.UseFreqList = 0;
			% Enable fmin, fmax, fstep
			enable_ui(handles.FminCtrl);
			enable_ui(handles.FmaxCtrl);
			enable_ui(handles.FstepCtrl);
			% set freqs, nfreqs to empty
			freqs = [];
			nfreqs = [];
		else
			%------------------------------------------------------
			% User selected file, so load freq list
			%------------------------------------------------------
			freqfile = fullfile(pathname, filename);
			[freqs, nfreqs] = NICal_LoadFreqList(freqfile);
		end
		
		if isempty(freqs)
			%------------------------------------------------------
			% user cancelled or freqs is empty for some reason, abort
			%------------------------------------------------------
			fprintf('%s: user cancelled\n', mfilename)
			% uncheck box, set UseFreqList in handles to 0
			update_ui_val(hObject, 0);
			handles.cal.UseFreqList = 0;
			% Enable fmin, fmax, fstep
			enable_ui(handles.FminCtrl);
			enable_ui(handles.FmaxCtrl);
			enable_ui(handles.FstepCtrl);
		else
			%------------------------------------------------------
			% valid freqs loaded, use them!
			%------------------------------------------------------
			% set UseFreqList to 1
			handles.cal.UseFreqList = 1;
			% disable fmin, fmax, fstep
			disable_ui(handles.FminCtrl);
			disable_ui(handles.FmaxCtrl);
			disable_ui(handles.FstepCtrl);
			% update Freqs, F and Nfreqs from freq file
			handles.cal.Freqs = freqs';
			handles.cal.Nfreqs = nfreqs;
			fprintf('\n\n*************\nFreqs:\n');
			fprintf('\t%.2f\n', handles.cal.Freqs);
			fprintf('*************\n\n');
		end
		
	else
		%------------------------------------------------------
		% user unchecked the box
		%------------------------------------------------------
		% set UseFreqList to 0
		handles.cal.UseFreqList = 0;
		% Enable fmin, fmax, fstep
		enable_ui(handles.FminCtrl);
		enable_ui(handles.FmaxCtrl);
		enable_ui(handles.FstepCtrl);
		% update Freqs and Nfreqs
		handles.cal.Freqs = handles.cal.Fmin:handles.cal.Fstep:handles.cal.Fmax;
		handles.cal.Nfreqs = length(handles.cal.Freqs);
	end
	guidata(hObject, handles)
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function FminCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	% check to make sure Fmin is in range [0.001 Fmax]
	if ~between(tmp, 1e-3, handles.cal.Fmax)
		% if not, revert to previous value and warn user with dialog box
		tmpstr = sprintf(['Min Freq must be greater than 0 ' ...
									'& less than Fmax (%d)'], handles.cal.Fmax);
		warndlg(	tmpstr, 'Invalid Min Freq');
		update_ui_str(hObject, handles.cal.Fmin);
	% check that Fmin + Fstep is in range
	elseif (tmp + handles.cal.Fstep) > handles.cal.Fmax
		tmpstr = 'Min Freq + Freq Step must be less than or equal to Fmax!';
		warndlg(	tmpstr, 'Invalid Min Freq');
		update_ui_str(hObject, handles.cal.Fmin);		
	% store the new Fmin value
	else
		handles.cal.Fmin = tmp;
		% update Freqs and Nfreqs
		handles.cal.Freqs = handles.cal.Fmin:handles.cal.Fstep:handles.cal.Fmax;
		handles.cal.Nfreqs = length(handles.cal.Freqs);
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function FmaxCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	% check range of Fmax [Fmin Fs/2]
	if ~between(tmp, handles.cal.Fmin,(handles.cal.Fs / 2))
		tmpstr = sprintf('Max Freq must be between Fmin (%d) & Fs / 2 (%f) Hz', ...
								handles.cal.Fmin, handles.cal.Fs);
		warndlg(tmpstr, 'Invalid Max Freq');
		update_ui_str(hObject, handles.cal.Fmax);
	% check that Fmin + Fstep range is still valid
	elseif (handles.cal.Fmin + handles.cal.Fstep) > tmp
		tmpstr = 'Min Freq + Freq Step must be less than or equal to Fmax!';
		warndlg(	tmpstr, 'Invalid Max Freq');
		update_ui_str(hObject, handles.cal.Fmax);		
	else
		% keep the new value
		handles.cal.Fmax = tmp;
		% update Freqs and Nfreqs
		handles.cal.Freqs = handles.cal.Fmin:handles.cal.Fstep:handles.cal.Fmax;
		handles.cal.Nfreqs = length(handles.cal.Freqs);
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function FstepCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	% step size needs to be in range [0.001 Fmax-Fmin]
	maxstep = handles.cal.Fmax - handles.cal.Fmin;
	if ~between(tmp, 0.001, maxstep)
		w = sprintf('FreqStep must be between 0.001 & Fmax-Fmin (%f)', maxstep);
		warndlg(w, 'Invalid FreqStep');
		update_ui_str(hObject, handles.cal.Fstep);
	else
		% store the new Fstep value
		handles.cal.Fstep = tmp;
		% update Freqs and Nfreqs
		handles.cal.Freqs = handles.cal.Fmin:handles.cal.Fstep:handles.cal.Fmax;
		handles.cal.Nfreqs = length(handles.cal.Freqs);
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function NrepsCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 1, 30)
		warndlg('# reps must be between 1 & 30', 'Invalid # Reps');
		update_ui_str(hObject, handles.cal.Nreps);
	else
		handles.cal.Nreps = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
% runs when AttenFixCtrl is checked or unchecked
%--------------------------------------------------------------------------
function AttenFixCtrl_Callback(hObject, eventdata, handles)
	%--------------------------------------------------------
	% set cal.AttenFix to the value of the control object
	%--------------------------------------------------------
	handles.cal.AttenFix = read_ui_val(hObject);
	%--------------------------------------------------------
	% update user interface depending on AttenFix value
	%--------------------------------------------------------
	if handles.cal.AttenFix
		%--------------------------------------------------------
		% use fixed attenuation 
		%--------------------------------------------------------
		% first, read the fixed attenuation value from AttenFixValueCtrl, store
		% in cal.AttenFixValue parameter
		handles.cal.AttenFixValue = read_ui_str(handles.AttenFixValueCtrl, 'n');
		% enable the AttenFixValueCtrl
		show_uictrl(handles.AttenFixValueCtrl);
		% make the text label for AttenFixValueCtrl visable
		set(handles.AttenFixValueCtrlText, 'Visible', 'on');
		% disable Minlevel, Maxlevel, Attenstep controls
		disable_ui(handles.MinlevelCtrl);
		disable_ui(handles.MaxlevelCtrl);
		disable_ui(handles.AttenStepCtrl);
		disable_ui(handles.StartAttenCtrl);
	else
		%--------------------------------------------------------
		% update settings for normal routine
		%--------------------------------------------------------
		% disable the AttenFixValueCtrl
		hide_uictrl(handles.AttenFixValueCtrl);
		% make the text label for AttenFixValueCtrl hidden
		set(handles.AttenFixValueCtrlText, 'Visible', 'off');
		% enable Minlevel, Maxlevel, Attenstep controls
		enable_ui(handles.MinlevelCtrl);
		enable_ui(handles.MaxlevelCtrl);
		enable_ui(handles.AttenStepCtrl);
		enable_ui(handles.StartAttenCtrl);
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function AttenFixValueCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(handles.AttenFixValueCtrl, 'n');
	if ~between(tmp, 0, 120)
		warndlg('fixed attenuation must be between 0 and 120', ...
						'Invalid Fixed Attenuation');
		update_ui_str(hObject, handles.cal.AttenFixValue);
	else
		handles.cal.AttenFixValue = tmp;
		guidata(hObject, handles);
	end	
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------
function MinlevelCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 0, handles.cal.Maxlevel)
		s = sprintf('Min Level must be between 0 & Max Level (%.1f) dB', ...
								handles.cal.Maxlevel);
		warndlg(s, 'Invalid Min Level');
		update_ui_str(hObject, handles.cal.Minlevel);
	else
		% store new Minlevel
		handles.cal.Minlevel = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function MaxlevelCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, handles.cal.Minlevel, 120)
		s = sprintf('Max Level must be between Min Level (%.1f) & 120 dB', ...
							handles.cal.Minlevel);
		warndlg(s, 'Invalid Max Level');
		update_ui_str(hObject, handles.cal.Maxlevel);
	else
		handles.cal.Maxlevel = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function AttenStepCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	maxstep = handles.cal.Maxlevel - handles.cal.Minlevel;
	if ~between(tmp, 0.1, maxstep)
		s = sprintf(['Atten Step must be between 0.1 & Max Level - ' ...
							'Min Level (%.1f) dB'], ...
							handles.cal.Maxlevel - handles.cal.Minlevel);
		warndlg(s, 'Invalid AttenStep');
		% reset Attenstep to original value
		update_ui_str(hObject, handles.cal.AttenStep);
	else
		% store new AttenStep value
		handles.cal.AttenStep = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function StartAttenCtrl_Callback(hObject, eventdata, handles)
	NICal_Constants;
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 0.1, MAX_ATTEN)
		s = sprintf('Start Atten must be between 0.1 & Max Atten (%.1f) dB', ...
							MAX_ATTEN);
		warndlg(s, 'Invalid StartAtten');
		% reset Attenstep to original value
		update_ui_str(hObject, handles.cal.StartAtten);
	else
		% store new AttenStep value
		handles.cal.StartAtten = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
function CollectBackgroundCtrl_Callback(hObject, eventdata, handles)
	handles.cal.CollectBackground = read_ui_val(hObject);
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% --- Executes on button press in MeasureLeakCtrl. ("Measure Crosstalk")
%--------------------------------------------------------------------------
function MeasureLeakCtrl_Callback(hObject, eventdata, handles)
	handles.cal.MeasureLeak = read_ui_val(hObject);
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% --- Executes on button press in SaveRawDataCtrl.
%--------------------------------------------------------------------------
function SaveRawDataCtrl_Callback(hObject, eventdata, handles)
	handles.cal.SaveRawData = read_ui_val(hObject);
	guidata(hObject, handles)
%--------------------------------------------------------------------------
	
%--------------------------------------------------------------------------
function CheckCalCtrl_Callback(hObject, eventdata, handles)
	handles.cal.CheckCal = read_ui_val(hObject)-1;
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% INPUT/OUTPUT SETTINGS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
function StimDurationCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 1, 1000)
		warndlg('Stimuli must be between 1 & 1000 ms', 'Invalid StimDuration');
		update_ui_str(hObject, handles.cal.StimDuration);
	else
		handles.cal.StimDuration = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function SweepDurationCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 10, 5000)
		warndlg('Sweep must be between 10 and 5000 ms', 'Invalid SweepDuration');
		update_ui_str(hObject, handles.cal.SweepDuration);
	else
		handles.cal.SweepDuration = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function StimDelayCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 0, 500)
		warndlg('StimDelay must be between 0 and 500 ms', 'Invalid StimDelay');
		update_ui_str(hObject, handles.cal.StimDelay);
	else
		handles.cal.StimDelay = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function StimRampCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, .1, 100)
		warndlg('StimRamp must be between .1 and 100 ms', 'Invalid StimRamp');
		update_ui_str(hObject, handles.cal.StimRamp);
	else
		handles.cal.StimRamp = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function ISICtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 0, 1000)
		warndlg('Interval between stimuli must be between 0 & 1000', ...
						'Invalid Interval');
		update_ui_str(hObject, handles.cal.ISI);
	else
		handles.cal.ISI = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% sets peak amplitude of output signal (tone)
%--------------------------------------------------------------------------
function DAscaleCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	% some checks first
	if isempty(newVal)
		warning('NICal:ValueOutOfRange', ...
					'%s: invalid tone level value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.DAscale);
	elseif ~all(isnumeric(newVal))
		warning('NICal:ValueOutOfRange',	...
						'%s: invalid tone level value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.DAscale);
	elseif ~between(newVal, 0, 10)
		warning('NICal:ValueOutOfRange', ...
							'%s: invalid tone level value %f', mfilename, newVal);
		disp('tone level must be between 0 and 10 Volts!')
		update_ui_str(hObject, handles.cal.DAscale);
	elseif length(newVal) ~= 1
		handles.cal.MicSensitivity = newVal(1);
		update_ui_str(hObject, handles.cal.DAscale);
	else
		% update cal.DAscale value to new value
		handles.cal.DAscale = newVal;
	end
	% update handles
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
% --- Executes on button press in InputFilterCtrl.
%-------------------------------------------------------------------------
function InputFilterCtrl_Callback(hObject, eventdata, handles)
	handles.cal.InputFilter = read_ui_val(hObject);
	if handles.cal.InputFilter
		% make filter freqs enabled
		show_uictrl(handles.HiPassFcCtrl);
		show_uictrl(handles.LoPassFcCtrl);
		set(handles.HiPassFcText, 'Visible', 'on');
		set(handles.LoPassFcText, 'Visible', 'on');
		set(handles.HiPassFcCtrlText, 'Visible', 'on');
		set(handles.LoPassFcCtrlText, 'Visible', 'on');
		% read values (might be redundant, but be sure)
		handles.cal.InputHPFc = read_ui_str(handles.HiPassFcCtrl, 'n');
		handles.cal.InputLPFc = read_ui_str(handles.LoPassFcCtrl, 'n');
	else
		hide_uictrl(handles.HiPassFcCtrl);
		hide_uictrl(handles.LoPassFcCtrl);
		set(handles.HiPassFcText, 'Visible', 'off');
		set(handles.LoPassFcText, 'Visible', 'off');
		set(handles.HiPassFcCtrlText, 'Visible', 'off');
		set(handles.LoPassFcCtrlText, 'Visible', 'off');
	end
		
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% sets bandpass filter lower Fc
%-------------------------------------------------------------------------
function HiPassFcCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	% some checks first
	if isempty(newVal)
		warndlg(sprintf('%s: invalid HiPass Fc %.2f', mfilename, newVal));
		update_ui_str(hObject, handles.cal.InputHPFc);
	elseif ~all(isnumeric(newVal))
		warning('NICal:ValueOutOfRange', ...
					'%s: invalid HiPassFc %.2f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.InputHPFc);
	elseif newVal <= 0
		warning('NICal:ValueOutOfRange', ...
					'%s: HiPassFc (%.2f) must be greater than 0', mfilename, newVal);
		update_ui_str(hObject, handles.cal.InputHPFc);
	elseif ~(newVal < handles.cal.InputLPFc)
		warning('NICal:ValueOutOfRange', ...
					'%s: HiPassFc (%.2f) must be less than LoPass Fc (%.2f)', ...
												mfilename, newVal, handles.cal.InputLPFc);
		update_ui_str(hObject, handles.cal.InputHPFc);
	elseif length(newVal) ~= 1
		handles.cal.InputHPFc = newVal(1);
		update_ui_str(hObject, handles.cal.InputHPFc);
	elseif newVal > (handles.cal.Fs / 2)
		warning('NICal:ValueOutOfRange', ...
					'%s: HiPassFc (%.2f) must be less than Fs/2 (%.2f)', ...
													mfilename, newVal, handles.cal.Fs);
		update_ui_str(hObject, handles.cal.InputHPFc);
	else
		% update cal.DAscale value to new value
		handles.cal.InputHPFc = newVal;
	end
	% update handles
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% runs when LoPassFc engaged
%-------------------------------------------------------------------------
function LoPassFcCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	%-------------------------------------------------
	% perform some checks on values
	%-------------------------------------------------
	if isempty(newVal)
		% newVal is empty
		warndlg(sprintf('%s: invalid LoPass Fc %.2f', mfilename, newVal));
		update_ui_str(hObject, handles.cal.InputLPFc);
	elseif ~all(isnumeric(newVal))
		% newval is not numeric
		warning('NICal:ValueOutOfRange', ...
						'%s: invalid Lo Pass Fc %.2f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.InputLPFc);
	elseif newVal <= 0
		% newval is less than or eq to zero
		warning('NICal:ValueOutOfRange', ...
					'%s: LoPassFc (%.2f) must be greater than 0', mfilename, newVal);
		update_ui_str(hObject, handles.cal.InputLPFc);
	elseif ~(newVal > handles.cal.InputHPFc)
		% newval is greater than hi pass value
		warning('NICal:ValueOutOfRange', ...
					'%s: LoPassFc (%.2f) must be greater than HiPass Fc (%.2f)', ...
												mfilename, newVal, handles.cal.InputHPFc);
		update_ui_str(hObject, handles.cal.InputLPFc);
	elseif newVal > (handles.cal.Fs / 2)
		warning('NICal:ValueOutOfRange', ...
					'%s: LoPassFc (%.2f) must be less than Fs/2 (%.2f)', ...
													mfilename, newVal, handles.cal.Fs);
		update_ui_str(hObject, handles.cal.InputLPFc);
	elseif length(newVal) ~= 1
		handles.cal.InputLPFc = newVal(1);
		update_ui_str(hObject, handles.cal.InputLPFc);
	else
		% update cal.DAscale value to new value
		handles.cal.InputLPFc = newVal;
	end
	% update handles
	guidata(hObject, handles);
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% MICROPHONE SETTINGS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Enables FR file use on button press in FRenableCtrl.
%--------------------------------------------------------------------------
function FRenableCtrl_Callback(hObject, eventdata, handles)
	currentState = read_ui_val(hObject);
	if currentState
		% do whatever needs to be done to enable mic FR
		enable_ui(handles.MicFRFileCtrl);
		disable_ui(handles.MicSensitivityCtrl);
		disable_ui(handles.DAscaleCtrl);
		update_ui_str(handles.MicFRFileCtrl, handles.cal.mic_fr_file);
		handles.cal.FRenable = 1;
	else
		% do whatever needs to be done to disable mic FR use
		disable_ui(handles.MicFRFileCtrl);
		enable_ui(handles.InputChannelCtrl);
		enable_ui(handles.MicGainCtrl);
		enable_ui(handles.MicSensitivityCtrl);
		enable_ui(handles.DAscaleCtrl);
		handles.cal.FRenable = 0;
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sets input channel on selection change in InputChannelCtrl.
%--------------------------------------------------------------------------
function InputChannelCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_val(hObject);
	if newVal == 1
		handles.cal.InputChannel = 1;
	elseif newVal == 2
		handles.cal.InputChannel = 2;
	elseif newVal == 3
		handles.cal.InputChannel = 3;
	else
		warning('NICal:ValueOutOfRange', '%s: unknown input channel %d', mfilename, newVal);
		update_ui_val(hObject, handles.cal.InputChannel);
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------
	
%--------------------------------------------------------------------------
% sets mic Gain
%--------------------------------------------------------------------------
function MicGainCtrl_Callback(hObject, eventdata, handles)
	% normally, read_ui_str(x, 'n') would be used, but str2double only
	% works on scalar inputs.  with 2 possible input channels, need to account
	% for both
	newVal = read_ui_str(hObject);
	newVal = str2num(newVal); %#ok<ST2NM>
	if isempty(newVal)
		warning('NICal:ValueOutOfRange', '%s: invalid MicGain value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.MicGain);
	elseif ~all(isnumeric(newVal))
		warning('NICal:ValueOutOfRange', '%s: invalid MicGain value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.MicGain);
% 	elseif length(newVal) ~= 1
% 		handles.cal.MicGain = newVal(1);
% 		update_ui_str(hObject, handles.cal.MicGain);
	else
		handles.cal.MicGain = newVal;
		% fix # gain values if Nchannels doesn't match # of gain values
		if handles.cal.Nchannels  ~= length(handles.cal.MicGain)
			handles.cal.MicGain = handles.cal.MicGain(1) .* ones(1, handles.cal.Nchannels);
			update_ui_str(handles.MicGainCtrl, handles.cal.MicGain);
		end
		% convert dB to linear scale
		handles.cal.Gain = invdb(handles.cal.MicGain);
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% sets Mic sensitiviy
%--------------------------------------------------------------------------
function MicSensitivityCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	if isempty(newVal)
		warning('NICal:ValueOutOfRange', '%s: invalid MicSensitivity value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.MicSensitivity);
	elseif ~all(isnumeric(newVal))
		warning('NICal:ValueOutOfRange', '%s: invalid MicSensitivity value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.MicSensitivity);
	elseif newVal <= 0
		warning('NICal:ValueOutOfRange', '%s: invalid MicSensitivity value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.MicSensitivity);
	elseif length(newVal) ~= 1
		handles.cal.MicSensitivity = newVal(1);
		update_ui_str(hObject, handles.cal.MicSensitivity);
	else
		handles.cal.MicSensitivity = newVal;
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% File settings control callbacks
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
% get microphone fr file for frdata
%--------------------------------------------------------------------------
function MicFRFileCtrl_Callback(hObject, eventdata, handles)

	tmpfile = read_ui_str(hObject);
	
	if ~exist(tmpfile, 'file')
		% user provided file does not exist
		h = warndlg('Microphone calibration file not found!', 'NICal Warning', 'modal');
		uiwait(h);
		
	
		[filename, pathname] = uigetfile( {'*.fr'; '*_fr.mat'}, ...
										'Load microphone calibration data from file...');
		if isequal(filename, 0) || isequal(pathname, 0)
			% user hit cancel, so restore the stored value of mic_fr_file
			update_ui_str(hObject, handles.cal.mic_fr_file);
			return
		else
			% set the tmp file to what user chose
			tmpfile = fullfile(pathname, filename);
		end
	end
	
	% now load fr data
	load(tmpfile, '-MAT', 'frdata');
	% make sure it was loaded
	if ~exist('frdata', 'var')
		warndlg('No frdata struct in fr file!', 'NICal Warning', 'modal');
		update_ui_str(hObject, handles.cal.mic_fr_file);
		guidata(hObject, handles);
		return
	else
		% store frdata and new file name
		handles.cal.mic_fr_file = tmpfile;
		handles.cal.mic_fr = frdata;
		guidata(hObject, handles);
		update_ui_str(handles.MicFRFileCtrl, handles.cal.mic_fr_file);
	end

%{

	% get the fr file data
	tmpfile = read_ui_str(handles.MicFRFileCtrl);
	if ~exist(tmpfile, 'file')
		h = warndlg('Microphone calibration file not found!', 'NICal Warning', 'modal');

		% revert to old value
		update_ui_str(handles.MicFRFileCtrl, handles.cal.mic_fr_file);
	else
		handles.cal.mic_fr_file = tmpfile;
		load(handles.cal.mic_fr_file, 'frdata');
		handles.cal.mic_fr = frdata;
		guidata(hObject, handles);
		update_ui_str(handles.MicFRFileCtrl, handles.cal.mic_fr_file);		
	end
%}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% get file to save data to
%--------------------------------------------------------------------------
function CalFileCtrl_Callback(hObject, eventdata, handles)
	oldfile = read_ui_str(hObject);
	[filename, pathname] = uiputfile('*.cal','Save calibration data to file', oldfile);
	if isequal(filename, 0) || isequal(pathname, 0)
		update_ui_str(hObject, handles.cal.calfile);
		return
	else
		handles.cal.calfile = fullfile(pathname, filename);
		update_ui_str(hObject, handles.cal.calfile);
		guidata(hObject, handles);
	end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over CalFileCtrl.
function CalFileCtrl_ButtonDownFcn(hObject, eventdata, handles)
	oldfile = read_ui_str(hObject);
	[filename, pathname] = uiputfile('*.cal','Save calibration data to file', oldfile);
	if isequal(filename, 0) || isequal(pathname, 0)
		update_ui_str(hObject, handles.cal.calfile);
		return
	else
		handles.cal.calfile = fullfile(pathname, filename);
		update_ui_str(hObject, handles.cal.calfile);
		guidata(hObject, handles);
	end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
function AutoSaveCtrl_Callback(hObject, eventdata, handles)
	handles.cal.AutoSave = read_ui_val(hObject);
	guidata(hObject, handles);
%--------------------------------------------------------------------------


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Menu Functions
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% File Menu
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function Menu_SaveCal_Callback(hObject, eventdata, handles)
	[calfile, calpath] = uiputfile('*.cal', ...
												'Save headphone calibration data in file');
	if calfile ~= 0
		% save the sequence so we can match up with the RF data
		datafile = fullfile(calpath, calfile);
		caldata = handles.caldata; %#ok<NASGU>
		save(datafile, '-MAT', 'caldata');
	end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function Menu_Close_Callback(hObject, eventdata, handles)
 	CloseRequestFcn(handles.figure1, eventdata, handles);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% MICROPHONE Menu
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function Menu_LoadFRData_Callback(hObject, eventdata, handles)
	% get the fr file data
	[frfile, frpath] = uigetfile('*.fr', 'Load FR data for earphones...');
	if frfile ~= 0
		handles.cal.mic_fr_file = fullfile(frpath, frfile);
		load(handles.cal.mic_fr_file, 'frdata', '-MAT');
		handles.cal.mic_fr = frdata;
		guidata(hObject, handles);
		update_ui_str(handles.MicFRFileCtrl, handles.cal.mic_fr_file);
	end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function Menu_CalibrateMicrophone_Callback(hObject, eventdata, handles)
	%---------------------------------------------------------------
	% turn off calibration ctrl, enable abort ctrl
	%---------------------------------------------------------------
	disable_ui(handles.RunCalibrationCtrl);
	show_uictrl(handles.AbortCtrl);
	set(handles.AbortCtrl, 'Value', 0);
	%---------------------------------------------------------------
	% initialize complete flag
	%---------------------------------------------------------------
	handles.CalComplete = 0;
	COMPLETE = 0; %#ok<NASGU>
	guidata(hObject, handles);
	%---------------------------------------------------------------
	% run microphone calibration script
 	%---------------------------------------------------------------
	NICal_MicrophoneCalibrate;
	%---------------------------------------------------------------
	% enable Calibration ctrl, disable abort ctrl
	%---------------------------------------------------------------
	enable_ui(handles.RunCalibrationCtrl);
	hide_uictrl(handles.AbortCtrl);
	set(handles.AbortCtrl, 'Value', 0);
	%---------------------------------------------------------------
	% save handles
	%---------------------------------------------------------------
	guidata(hObject, handles);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Analysis Menu
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function Menu_ProcessTriggeredData_Callback(hObject, eventdata, handles)
	disp('This will process .dat data collected during triggered acquisition')
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Calibration Settings Menu
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function Menu_LoadSettings_Callback(hObject, eventdata, handles)
	% get the settings file name
	[sfilename, sfilepath] = uigetfile(	'*_settings.mat', ...
													'Load Calibration Settings...', ...
													handles.userconfigpath );
	if sfilename ~= 0
		cal = [];
		load(fullfile(sfilepath, sfilename), 'cal');
		handles.cal = cal;
		guidata(hObject, handles);
	end
	% update user interface
	NICal_UpdateUIFromCal(handles, handles.cal);		
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function Menu_SaveSettings_Callback(hObject, eventdata, handles)
	% get the settings file name
	[sfilename, sfilepath] = uiputfile(	'*_settings.mat', ...
													'Save Calibration Settings...', ...
													handles.userconfigpath );
	if sfilename ~= 0
		cal = handles.cal; %#ok<NASGU>
		save(fullfile(sfilepath, sfilename), '-MAT', 'cal');
	end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function Menu_EditCalibrationSettings_Callback(hObject, eventdata, handles)
	% make local copy of cal
	existing_cal = handles.cal;
	% open StructDlg to edit cal
	try
		new_cal = StructDlg(	existing_cal, ...
									'Cal Struct Settings', ...
									[], [], [], [], ...
									'vert_spacing', 0.5, ...
									'font_size', 8 );
	catch errMsg
		save('NICal.err', 'errMsg', '-MAT');
		errMsg.stack
		new_cal = [];
	end
	
	if isempty(new_cal)
		% user closed window so do nothing
		disp('Reverting to original cal settings...')
		return
	else
		handles.cal = new_cal;
		guidata(hObject, handles);
		% update user interface
		NICal_UpdateUIFromCal(handles, handles.cal);
	end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Load default cal settings
%--------------------------------------------------------------------------
function Menu_ReloadDefaults_Callback(hObject, eventdata, handles)
	fprintf('Reloading cal settings from defaults file %s ...\n', ...
														handles.defaultsfile)
	load(handles.defaultsfile, 'cal');
	handles.cal = cal;
	clear cal
	guidata(hObject, handles);
	% update user interface
	NICal_UpdateUIFromCal(handles, handles.cal);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% save cal settings as default
%--------------------------------------------------------------------------
function Menu_SaveAsDefaultSettings_Callback(hObject, eventdata, handles)
	s = sprintf('Saving cal defaults in file %s', handles.defaultsfile);
	msgbox(s, 'NICal: Save Defaults', 'non-modal');
	cal = handles.cal; %#ok<NASGU>
	uical = NICal_UpdateCalFromUI(handles); %#ok<NASGU>
	save('calAB.mat', 'cal', 'uical', '-MAT');
	save(handles.defaultsfile, 'cal');
	clear cal
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% DEBUG Menu
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function Menu_DumpHandles_Callback(hObject, eventdata, handles)
	save('NICalhandles.mat', 'handles', '-MAT')
%--------------------------------------------------------------------------
function Menu_DumpNISettings_Callback(hObject, eventdata, handles)
	% Load the settings and constants
	NICal_Constants;
	NICal_settings;
	% make a local copy of the cal settings structure
	cal = handles.cal; %#ok<NASGU>
	% save the GUI handle information
	guidata(hObject, handles);
	% make local copy of iodev struct
	iodev = handles.iodev; %#ok<NASGU>
	handles.Nchannels = 2;
	guidata(hObject, handles);
	% Start DAQ things
	NICal_NIinit;
	guidata(hObject, handles);
	% save handles
	fname = fullfile(pwd, 'aiao.mat');
	[fname, pname] = uiputfile('*.mat', 'Save ai and ao structs to file', fname);
	if isequal(fname, 0) || isequal(pname, 0)
		fname = 'aiao.mat';
		pname = pwd;
	end
	ai = handles.iodev.NI.ai; %#ok<NASGU>
	ao = handles.iodev.NI.ao; %#ok<NASGU>
	iodev = handles.iodev; %#ok<NASGU>
	save(fullfile(pname, fname), 'ai', 'ao', 'iodev', '-MAT');
	clear ai ao
	% Stop DAQ
	NICal_NIexit;
	guidata(hObject, handles);
%--------------------------------------------------------------------------
function Menu_Keyboard_Callback(hObject, eventdata, handles)
	keyboard
%-------------------------------------------------------------------------
function Menu_DAQReset_Callback(hObject, eventdata, handles)
	if handles.DAQSESSION
		try
			daqreset
		catch errMsg
			fprintf('\nProblem with daqreset!\n\n')
			disp(errMsg)
			return
		end
	else
		warning('%s: not using NIDAQ Session interface', mfilename);
	end			
%-------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Options Menu
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function Menu_ToneStack_Callback(hObject, eventdata, handles)
	newVal = get(handles.Menu_ToneStack, 'Checked');
	if strcmpi(newVal, 'off')
		set(handles.Menu_ToneStack, 'Checked', 'on');
		handles.ToneStack = 1;
		set(handles.Menu_ToneSweep, 'Checked', 'off');
		handles.ToneSweep = 0;
	else
		set(handles.Menu_ToneStack, 'Checked', 'off');
		handles.ToneStack = 0;
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------
function Menu_ToneSweep_Callback(hObject, eventdata, handles)
	newVal = get(handles.Menu_ToneSweep, 'Checked');
	if strcmpi(newVal, 'off')
		set(handles.Menu_ToneSweep, 'Checked', 'on');
		handles.ToneSweep = 1;
		set(handles.Menu_ToneStack, 'Checked', 'off');
		handles.ToneStack = 0;
	else
		set(handles.Menu_ToneSweep, 'Checked', 'off');
		handles.ToneStack = 0;
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------
function Menu_ContinuousRecord_Callback(hObject, eventdata, handles)
	newVal = get(handles.Menu_ContinuousRecord, 'Checked');
	if strcmpi(newVal, 'off')
		set(handles.Menu_ToneSweep, 'Checked', 'off');
		handles.ToneSweep = 0;
		set(handles.Menu_ToneStack, 'Checked', 'off');
		handles.ToneStack = 0;
		set(handles.Menu_ContinuousRecord, 'Checked', 'on');
		handles.ContinuousRecord = 1;
	else
		set(handles.Menu_ContinuousRecord, 'Checked', 'off');
		handles.ContinuousRecord = 0;
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% GUI I/O, misc Functions
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = NICal_OutputFcn(hObject, eventdata, handles) 
	% Get default command line output from handles structure
	varargout{1} = [];
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function CloseRequestFcn(hObject, eventdata, handles)
	pause(0.1);
	delete(hObject);
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Create Functions
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function bgcolorset(hObject)
	if ispc && isequal(get(hObject,'BackgroundColor'), ...
												get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
%--------------------------------------------------------------------------
function FreqValText_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
	bgcolorset(hObject);
function LValText_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function LSPLText_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function RValText_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function RSPLText_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function FminCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function FmaxCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function FstepCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function MinlevelCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function MaxlevelCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function AttenStepCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function NrepsCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function ISICtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function LAttenText_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function RAttenText_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function StartAttenCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function SideCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function CheckCalCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function AttenFixValueCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function MicFRFileCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function InputChannelCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function MicGainCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function MicSensitivityCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function DAscaleCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function CalFileCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function HiPassFcCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function LoPassFcCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function StimDurationCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function SweepDurationCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function StimDelayCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
function StimRampCtrl_CreateFcn(hObject, eventdata, handles)
	bgcolorset(hObject);
%-------------------------------------------------------------------------



