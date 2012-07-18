function cal = UpdateCalFromUI(handles)
%--------------------------------------------------------------------------
% NICal_UpdateUIFromCal(handles, cal)
%--------------------------------------------------------------------------
% NICal program
% TytoLogy Project
%------------------------------------------------------------------------
% updates CAL from settings in UI controls
%------------------------------------------------------------------------
% Input Arguments:
%	handles				handles struct from GUI
%
% Output Arguments:
% 	cal					cal struct
%------------------------------------------------------------------------
% See also: NICal, NICal_UpdateUIFromCal, NICal_calstruct_init
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 18 July 2012 (SJS)
% 
% Revisions:
%--------------------------------------------------------------------------

%-----------------------------------------
% CALIBRATION SETTINGS
%-----------------------------------------
cal.Side = read_ui_val(handles.SideCtrl);
cal.Nreps = read_ui_str(handles.NrepsCtrl, 'n');
% Freq settings
cal.UseFreqList = read_ui_val(handles.FreqListCtrl);
cal.Fmin = read_ui_str(handles.FminCtrl, 'n');
cal.Fmax = read_ui_str(handles.FmaxCtrl, 'n');
cal.Fstep = read_ui_str(handles.FstepCtrl, 'n');
% update Freqs and Nfreqs
cal.Freqs = cal.Fmin:cal.Fstep:cal.Fmax;
cal.Nfreqs = length(cal.Freqs);
% Attenuation settings
cal.AttenFix = read_ui_val(handles.AttenFixCtrl);
cal.AttenFixValue = read_ui_str(handles.AttenFixValueCtrl, 'n');
cal.Minlevel = read_ui_str(handles.MinlevelCtrl, 'n');
cal.Maxlevel = read_ui_str(handles.MaxlevelCtrl, 'n');
cal.AttenStep = read_ui_str(handles.AttenStepCtrl, 'n');
cal.StartAtten = read_ui_str(handles.StartAttenCtrl, 'n');
%  Check Cal setting
cal.CheckCal = read_ui_val(handles.CheckCalCtrl) - 1;
%-----------------------------------------
% INPUT/OUTPUT SETTINGS
%-----------------------------------------
cal.ISI = read_ui_str(handles.ISICtrl, 'n');
cal.StimDuration = read_ui_str(handles.StimDurationCtrl, 'n');
cal.StimDelay = read_ui_str(handles.StimDelayCtrl, 'n');
cal.StimRamp = read_ui_str(handles.StimRampCtrl, 'n');
cal.SweepDuration = read_ui_str(handles.SweepDurationCtrl, 'n');
cal.DAscale = read_ui_str(handles.DAscaleCtrl, 'n');
cal.MeasureLeak = read_ui_val(handles.MeasureLeakCtrl);
% some other derived variables (not user-settable)
% Total time to acquire data (ms)
cal.AcqDuration = cal.SweepDuration;
% Total sweep time = sweep duration + inter stimulus interval (ms)
cal.SweepPeriod = cal.SweepDuration + cal.ISI;
% input bandpass filter
cal.InputFilter = read_ui_val(handles.InputFilterCtrl);
cal.InputHPFc = read_ui_str(handles.HiPassFcCtrl, 'n');
cal.InputLPFc = read_ui_str(handles.LoPassFcCtrl, 'n');
cal.Fs = handles.cal.Fs;
% Nyquist frequency
fnyq = cal.Fs/2;
% filter order
cal.forder = 5;
% passband definition
cal.fband = [cal.InputHPFc cal.InputLPFc] ./ fnyq;
% filter coefficients using a butterworth bandpass filter
[cal.fcoeffb, cal.fcoeffa] = butter(cal.forder, cal.fband, 'bandpass');
%-----------------------------------------
% MICROPHONE SETTINGS
%-----------------------------------------
cal.FRenable = read_ui_val(handles.FRenableCtrl);
cal.InputChannel = read_ui_val(handles.InputChannelCtrl);
cal.MicGain = read_ui_str(handles.MicGainCtrl, 'n');
cal.MicSensitivity = read_ui_str(handles.MicSensitivityCtrl, 'n');		
%-----------------------------------------
% FILE SETTINGS
%-----------------------------------------
cal.mic_fr_file = read_ui_str(handles.MicFRFileCtrl);
cal.calfile = read_ui_str(handles.CalFileCtrl);
cal.AutoSave = read_ui_val(handles.AutoSaveCtrl);
