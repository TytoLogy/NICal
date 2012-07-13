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

% Last Modified by GUIDE v2.5 13-Jul-2012 15:47:40

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
function NICal_OpeningFcn(hObject, eventdata, handles, varargin)
	%----------------------------------------------------------
	%----------------------------------------------------------
	% Setup Paths
	%----------------------------------------------------------
	%----------------------------------------------------------
	disp([mfilename ': checking paths'])
	pdir = ['C:\TytoLogy\TytoSettings\' getenv('USERNAME')];
	if isempty(which('RPload'))
		% could not find the RPload.m function (which is in TytoLogy
		% toolbox) which suggests that the paths are not set or are 
		% incorrect for this setup.  load the paths using the tytopaths program.
		%--------
		% First, store the current path
		cdir = pwd;
		% build the path to the user's TytoSettings directory and
		% change dirs to it.  Run the tytopaths script and then
		% return to the original ("current") directory
		disp([mfilename ': loading paths using ' pdir])
		cd(pdir);
		tytopaths
		cd(cdir);
	else
		disp([mfilename ': paths ok, launching programn'])
	end
	
	%----------------------------------------------------------
	%----------------------------------------------------------
	% load the configuration information, store in config structure
	% The HPSearch_Configuration.m function file will usually live in the
	% <tytology path>\TytoSettings\<username\ directory
	%----------------------------------------------------------
	%----------------------------------------------------------
	% load the configuration information, store in config structure
	if isempty(which('NICal_Configuration'))
		% need to add user config path
		addpath(['C:\TytoLogy\TytoSettings\' getenv('USERNAME')]);
	end	
	config = NICal_Configuration;
	% save handles
	guidata(hObject, handles);	
	
	%----------------------------------------------------------
	%----------------------------------------------------------
	% Initial Calibration settings
	%----------------------------------------------------------
	%----------------------------------------------------------
	
	%------------------------------------------------------------
	% first check to see if defaults file exists
	%------------------------------------------------------------
	handles.defaultsfile = fullfile(pdir, [mfilename '_Defaults.mat']);
	if exist(handles.defaultsfile, 'file')
		fprintf('Loading cal settings from defaults file %s ...\n', handles.defaultsfile)
		load(handles.defaultsfile, 'cal');
	else
		% no defaults found, so use internal values
		disp('no defaults found, using internal values')
		handles.cal = cal_structinit;
	end
	guidata(hObject, handles);

	%----------------------------------------------------------
	% update user interface
	%----------------------------------------------------------
	UpdateUIFromCal(handles, handles.cal);
	
	%----------------------------------------------------------
	%----------------------------------------------------------
	% setup iodev NI device interface struct
	% using information from the config struct 
	%----------------------------------------------------------
	%----------------------------------------------------------
	iodev.Circuit_Path = config.CIRCUIT_PATH;
	iodev.Circuit_Name = config.CIRCUIT_NAME;
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
	handles.initfunction = config.IOINITFUNCTION;
	handles.iofunction = config.IOFUNCTION;
	handles.attfunction = config.ATTENFUNCTION;
	guidata(hObject, handles);
	
	%----------------------------------------------------------
	%----------------------------------------------------------
	% setup for not using FR file
	%----------------------------------------------------------
	%----------------------------------------------------------
	disable_ui(handles.MicFRFileCtrl);
	%%% whatever needs to be done to disable mic FR use
	disable_ui(handles.MicFRFileCtrl);
	enable_ui(handles.InputChannelCtrl);
	enable_ui(handles.MicGainCtrl);
	enable_ui(handles.MicSensitivityCtrl);
	enable_ui(handles.DAscaleCtrl);
	handles.FRenable = 0;
	guidata(hObject, handles);
	
	%----------------------------------------------------------
	%----------------------------------------------------------
	% default is no leak (crosstalk) measurement
	%----------------------------------------------------------
	%----------------------------------------------------------
	handles.MeasureLeak = 0;
	update_ui_val(handles.MeasureLeakCtrl, handles.MeasureLeak);
	guidata(hObject, handles);
	
	%----------------------------------------------------------
	%----------------------------------------------------------
	% Update handles structure
	%----------------------------------------------------------
	%----------------------------------------------------------
	handles.CalComplete = 0;
	handles.output = hObject;
	guidata(hObject, handles);		
		
	%----------------------------------------------------------
	%----------------------------------------------------------
	% UIWAIT makes NICal wait for user response (see UIRESUME)
	%----------------------------------------------------------
	%----------------------------------------------------------
	% uiwait(handles.figure1);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% updates UI from settings in cal
%--------------------------------------------------------------------------
function UpdateUIFromCal(handles, cal)
	update_ui_str(handles.Fmin, cal.Fmin);
	update_ui_str(handles.Fmax, cal.Fmax);
	update_ui_str(handles.Fstep, cal.Fstep);
	update_ui_str(handles.Minlevel, cal.Minlevel);
	update_ui_str(handles.Maxlevel, cal.Maxlevel);
	update_ui_str(handles.AttenStep, cal.AttenStep);
	update_ui_str(handles.Nreps, cal.Nreps);
	update_ui_str(handles.ISI, cal.ISI);
	update_ui_val(handles.Side, cal.Side);
	update_ui_val(handles.AutoSave, cal.AutoSave);
	update_ui_val(handles.CheckCalCtrl, cal.CheckCal + 1);
	update_ui_val(handles.AttenFixCtrl, cal.AttenFix);
	update_ui_str(handles.AttenFixValueCtrl, cal.AttenFixValue);
	set(handles.AttenFixCtrl, 'HitTest', 'on');
	set(handles.AttenFixCtrl, 'Enable', 'on');
	set(handles.AttenFixCtrl, 'Visible', 'on');
	if cal.AttenFix
		set(handles.AttenFixValueCtrl, 'HitTest', 'on');
		set(handles.AttenFixValueCtrl, 'Enable', 'on');
		set(handles.AttenFixValueCtrl, 'Visible', 'on');
		set(handles.AttenFixValueCtrlText, 'HitTest', 'off');
		set(handles.AttenFixValueCtrlText, 'Enable', 'on');
		set(handles.AttenFixValueCtrlText, 'Visible', 'on');
	else
		set(handles.AttenFixValueCtrl, 'HitTest', 'off');
		set(handles.AttenFixValueCtrl, 'Enable', 'on');
		set(handles.AttenFixValueCtrl, 'Visible', 'off');
		set(handles.AttenFixValueCtrlText, 'HitTest', 'off');
		set(handles.AttenFixValueCtrlText, 'Enable', 'off');
		set(handles.AttenFixValueCtrlText, 'Visible', 'off');
	end
	update_ui_str(handles.MicFRFileCtrl, cal.mic_fr_file);

	update_ui_val(handles.InputChannelCtrl, cal.InputChannel);
	update_ui_str(handles.MicGainCtrl, cal.MicGain);
	update_ui_str(handles.MicSensitivityCtrl, cal.MicSensitivity);
	update_ui_str(handles.DAscaleCtrl, cal.DAscale);
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Run Calibration callback
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function RunCalibration_ctrl_Callback(hObject, eventdata, handles)
	%---------------------------------------------------------------
	% turn off calibration ctrl, enable abort ctrl
	%---------------------------------------------------------------
	set(handles.RunCalibration_ctrl, 'Enable', 'off');
	set(handles.Abort_ctrl, 'Enable', 'on');
	set(handles.Abort_ctrl, 'Visible', 'on');
	set(handles.Abort_ctrl, 'HitTest', 'on');
	set(handles.Abort_ctrl, 'Value', 0);
	
	%---------------------------------------------------------------
	% initialize complete flag
	%---------------------------------------------------------------
	handles.CalComplete = 0;
	COMPLETE = 0;
	guidata(hObject, handles);

	%---------------------------------------------------------------
	% run calibration script
 	%---------------------------------------------------------------
	NICal_RunCalibration

	%---------------------------------------------------------------
	% enable Calibration ctrl, disable abort ctrl
	%---------------------------------------------------------------
	set(handles.RunCalibration_ctrl, 'Enable', 'on');
	set(handles.Abort_ctrl, 'Enable', 'off');
	set(handles.Abort_ctrl, 'Visible', 'off');
	set(handles.Abort_ctrl, 'HitTest', 'off');
	set(handles.Abort_ctrl, 'Value', 0);
	
	%---------------------------------------------------------------
	% save cal settings
	%---------------------------------------------------------------
	if COMPLETE
		handles.CalComplete = 1;
		save(handles.defaultsfile, 'cal');
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Save Calibration Button
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Save the data 
function SaveCalibration_Callback(hObject, eventdata, handles)
	if handles.CalComplete
		[calfile, calpath] = uiputfile('*_cal.mat','Save headphone calibration data in file');
		if calfile ~= 0
			% save the sequence so we can match up with the RF data
			datafile = fullfile(calpath, calfile);
			caldata = handles.caldata;
			save(datafile, '-MAT', 'caldata');
		end
	end
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Calibration Settings callbacks
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%{
function FreqVal_Callback(hObject, eventdata, handles)
function LValText_Callback(hObject, eventdata, handles)
function LSPLText_Callback(hObject, eventdata, handles)
function RValText_Callback(hObject, eventdata, handles)
function RSPLText_Callback(hObject, eventdata, handles)
function LAttenText_Callback(hObject, eventdata, handles)
function RAttenText_Callback(hObject, eventdata, handles)
%}
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function Fmin_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 0, handles.cal.Fmax)
		tmpstr = sprintf('Min Freq must be between 0 & Fmax (%d)', handles.cal.Fmax);
		warndlg(	tmpstr, 'Invalid Min Freq');
		update_ui_str(hObject, handles.cal.Fmin);
	else
		handles.cal.Fmin = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
	
%-------------------------------------------------------------------------
function Fmax_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, handles.cal.Fmin, 125000)
		tmpstr = sprintf('Max Freq must be between Fmin (%d) & 125,000 Hz', ...
								handles.cal.Fmin);
		warndlg(tmpstr, 'Invalid Max Freq');
		update_ui_str(hObject, handles.cal.Fmax);
	else
		handles.cal.Fmax = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function Fstep_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	maxstep = handles.cal.Fmax - handles.cal.Fmin;
	if ~between(tmp, 1, maxstep)
		warndlg('FreqStep must be between 1 & Fmax-Fmin', 'Invalid FreqStep');
		update_ui_str(hObject, handles.cal.Fstep);
	else
		handles.cal.Fstep = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function Minlevel_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 1, handles.cal.Maxlevel)
		warndlg('Min Level must be between 0 & Max Level', 'Invalid Min Level');
		update_ui_str(hObject, handles.cal.Minlevel);
	else
		handles.cal.Minlevel = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function Maxlevel_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, handles.cal.Minlevel, 120)
		warndlg('Max Level must be between Min Level & 120 dB', 'Invalid Max Level');
		update_ui_str(hObject, handles.cal.Maxlevel);
	else
		handles.cal.Maxlevel = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function AttenStep_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	maxstep = handles.cal.Maxlevel - handles.cal.Minlevel;
	if ~between(tmp, 1, maxstep)
		warndlg('Atten Step must be between 1 & Max Level - Min Level', 'Invalid AttenStep');
		update_ui_str(hObject, handles.cal.AttenStep);
	else
		handles.cal.AttenStep = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------
	
%-------------------------------------------------------------------------
function Nreps_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 1, 30)
		warndlg('# reps must be between 1 & 30', 'Invalid # Reps');
		update_ui_str(hObject, handles.cal.Nreps);
	else
		handles.cal.Nreps = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function ISI_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(hObject, 'n');
	if ~between(tmp, 0, 1000)
		warndlg('Interval between stimuli must be between 0 & 1000', 'Invalid Interval');
		update_ui_str(hObject, handles.cal.ISI);
	else
		handles.cal.ISI = tmp;
		guidata(hObject, handles);
	end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
function AutoSave_Callback(hObject, eventdata, handles)
	handles.cal.AutoSave = read_ui_val(hObject);
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function Abort_ctrl_Callback(hObject, eventdata, handles)
	disp('ABORTING Calibration!')
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function CheckCalCtrl_Callback(hObject, eventdata, handles)
	handles.cal.CheckCal = read_ui_val(hObject)-1;
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function AttenFixValueCtrl_Callback(hObject, eventdata, handles)
	tmp = read_ui_str(handles.AttenFixValueCtrl, 'n');
	if ~between(tmp, 0, 120)
		warndlg('fixed attenuation must be between 0 and 120', 'Invalid Fixed Attenuation');
		update_ui_str(hObject, handles.cal.AttenFixValue);
	else
		handles.cal.AttenFixValue = tmp;
		guidata(hObject, handles);
	end	
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function AttenFixCtrl_Callback(hObject, eventdata, handles)
	handles.cal.AttenFix = read_ui_val(hObject);
	% update user interface
	if handles.cal.AttenFix
		% update settings for fixed attenuation routine
		set(handles.AttenFixValueCtrl, 'HitTest', 'on');
		set(handles.AttenFixValueCtrl, 'Enable', 'on');
		set(handles.AttenFixValueCtrl, 'Visible', 'on');
		set(handles.AttenFixValueCtrlText, 'HitTest', 'off');
		set(handles.AttenFixValueCtrlText, 'Enable', 'on');
		set(handles.AttenFixValueCtrlText, 'Visible', 'on');
		handles.cal.AttenFixValue = read_ui_str(handles.AttenFixValueCtrl, 'n');
	else
		% update settings for normal routine
		set(handles.AttenFixValueCtrl, 'HitTest', 'off');
		set(handles.AttenFixValueCtrl, 'Enable', 'on');
		set(handles.AttenFixValueCtrl, 'Visible', 'off');
		set(handles.AttenFixValueCtrlText, 'HitTest', 'off');
		set(handles.AttenFixValueCtrlText, 'Enable', 'off');
		set(handles.AttenFixValueCtrlText, 'Visible', 'off');
	end
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function MicFRFileCtrl_Callback(hObject, eventdata, handles)
	% get the fr file data
	tmpfile = read_ui_str(handles.MicFRFileCtrl);
	if ~exist(tmpfile, 'file')
		warndlg('Microphone calibration file not found!', 'NICal Warning');
		% revert to old value
		update_ui_str(handles.MicFRFileCtrl, handles.cal.mic_fr_file);
	else
		handles.cal.mic_fr_file = tmpfile;
		load(handles.cal.mic_fr_file, 'frdata');
		handles.cal.mic_fr = frdata;
		guidata(hObject, handles);
		update_ui_str(handles.MicFRFileCtrl, handles.cal.mic_fr_file);		
	end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function CalFileCtrl_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% --- Executes on button press in MeasureLeakCtrl.
%--------------------------------------------------------------------------
function MeasureLeakCtrl_Callback(hObject, eventdata, handles)
	handles.MeasureLeak = read_ui_val(hObject);
	guidata(hObject, handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Enables FR file use on button press in FRenableCtrl.
%--------------------------------------------------------------------------
function FRenableCtrl_Callback(hObject, eventdata, handles)
	currentState = read_ui_val(hObject);
	if currentState
		% do whatever needs to be done to enable mic FR
		enable_ui(handles.MicFRFileCtrl);
		disable_ui(handles.InputChannelCtrl);
		disable_ui(handles.MicGainCtrl);
		disable_ui(handles.MicSensitivityCtrl);
		disable_ui(handles.DAscaleCtrl);
		update_ui_str(handles.MicFRFileCtrl, handles.cal.mic_fr_file);
		handles.FRenable = 1;
	else
		disable_ui(handles.MicFRFileCtrl);
		% do whatever needs to be done to disable mic FR use
		disable_ui(handles.MicFRFileCtrl);
		enable_ui(handles.InputChannelCtrl);
		enable_ui(handles.MicGainCtrl);
		enable_ui(handles.MicSensitivityCtrl);
		enable_ui(handles.DAscaleCtrl);
		handles.FRenable = 0;
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
	else
		handles.cal.InputChannel = 2;
	end
%--------------------------------------------------------------------------
	
%--------------------------------------------------------------------------
% sets mic Gain
%--------------------------------------------------------------------------
function MicGainCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	
	if isempty(newVal)
		warning('NICal:ValueOutOfRange', '%s: invalid MicGain value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.MicGain);
	elseif ~all(isnumeric(newVal))
		warning('NICal:ValueOutOfRange', '%s: invalid MicGain value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.MicGain);
	elseif length(newVal) ~= 1
		handles.cal.MicGain = newVal(1);
		update_ui_str(hObject, handles.cal.MicGain);
	else
		handles.cal.MicGain = newVal;
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

%--------------------------------------------------------------------------
% sets peak amplitude of output signal (tone)
%--------------------------------------------------------------------------
function DAscaleCtrl_Callback(hObject, eventdata, handles)
	newVal = read_ui_str(hObject, 'n');
	% some checks first
	if isempty(newVal)
		warning('NICal:ValueOutOfRange', '%s: invalid tone level value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.DAscale);
	elseif ~all(isnumeric(newVal))
		warning('NICal:ValueOutOfRange', '%s: invalid  tone level value %f', mfilename, newVal);
		update_ui_str(hObject, handles.cal.DAscale);
	elseif ~between(newVal, 0, 10)
		warning('NICal:ValueOutOfRange', '%s: invalid  tone level value %f', mfilename, newVal);
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
%-------------------------------------------------------------------------



%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function HiPassFcCtrl_Callback(hObject, eventdata, handles)





%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function LoPassFcCtrl_Callback(hObject, eventdata, handles)



%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Menu Functions
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
function Menu_SaveCal_Callback(hObject, eventdata, handles)
	[calfile, calpath] = uiputfile('*_cal.mat','Save headphone calibration data in file');
	if calfile ~= 0
		% save the sequence so we can match up with the RF data
		datafile = fullfile(calpath, calfile);
		caldata = handles.caldata;
		save(datafile, '-MAT', 'caldata');
	end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function Menu_LoadFRData_Callback(hObject, eventdata, handles)
	% get the fr file data
	[frfile, frpath] = uigetfile('*_fr.mat', 'Load FR data for earphones...');
	if frfile ~= 0
		handles.cal.mic_fr_file = fullfile(frpath, frfile);
		load(handles.cal.mic_fr_file, 'frdata');
		handles.cal.mic_fr = frdata;
		guidata(hObject, handles);
		update_ui_str(handles.MicFRFileCtrl, handles.cal.mic_fr_file);
	end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function Menu_Close_Callback(hObject, eventdata, handles)
 	CloseRequestFcn(handles.figure1, eventdata, handles);
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
% TDT menu (unused)
%-------------------------------------------------------------------------
%{
function TDTSettingsMenuCtrl_Callback(hObject, eventdata, handles)
	iodev = handles.iodev;
	fullcircuit = fullfile(iodev.Circuit_Path, [iodev.Circuit_Name '.rcx'])
	if ~exist(fullcircuit, 'file')
		warning('%s: circuit %s not found...', mfilename, fullcircuit)
		[fname, pname] = uigetfile('*.rcx', 'Select TDT RPvD circuit file');
		if fname == 0
			% user cancelled request
			return
		end
		% need to strip off .rcx from filename
		[tmp1, fname, fext] = fileparts(fname)
		iodev.Circuit_Name = fname;
		iodev.Circuit_Path = pname;
		handles.iodev = iodev;
		guidata(hObject, handles);
		iodev
		
	else
		[fname, pname] = uigetfile('*.rcx', ...
								'Select TDT RPvD circuit file', ...
								fullcircuit);
		if fname == 0
			% user cancelled request
			return
		end
		% need to strip off .rcx from filename
		[tmp1, fname, fext] = fileparts(fname);
		iodev.Circuit_Name = fname;
		iodev.Circuit_Path = pname;
		handles.iodev = iodev;
		guidata(hObject, handles);
		iodev
	end
%}
%-------------------------------------------------------------------------



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
	varargout{1} = handles.output;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function CloseRequestFcn(hObject, eventdata, handles)
	pause(0.1);
	delete(hObject);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% --- Executes on selection change in Side.
function Side_Callback(hObject, eventdata, handles)
	newVal = read_ui_val(hObject);
	handles.cal.Side = newVal;
	guidata(hObject, handles);
%--------------------------------------------------------------------------
	



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Create Functions
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function FreqValText_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function LValText_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function LSPLText_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function RValText_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function RSPLText_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function Fmin_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function Fmax_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function Fstep_CreateFcn(hObject, eventdata, handles)
function Minlevel_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function Maxlevel_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function AttenStep_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function Nreps_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function ISI_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function LAttenText_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function RAttenText_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		set(hObject,'BackgroundColor','white');
	end
function Side_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CheckCalCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function AttenFixValueCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function MicFRFileCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	    set(hObject,'BackgroundColor','white');
	end
function InputChannelCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function MicGainCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function MicSensitivityCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function DAscaleCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function CalFileCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
function HiPassFcCtrl_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LoPassFcCtrl_CreateFcn(hObject, eventdata, handles)
	if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
		 set(hObject,'BackgroundColor','white');
	end
%-------------------------------------------------------------------------



