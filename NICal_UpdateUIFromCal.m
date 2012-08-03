function NICal_UpdateUIFromCal(handles, cal)
%--------------------------------------------------------------------------
% NICal_UpdateUIFromCal(handles, cal)
%--------------------------------------------------------------------------
% NICal program
% TytoLogy Project
%------------------------------------------------------------------------
% updates UI features from settings in cal structure
%------------------------------------------------------------------------
% Input Arguments:
%	handles				handles struct from GUI
% 	cal					cal struct
%
% Output Arguments:
%	-NONE-
%------------------------------------------------------------------------
% See also: NICal, NICal_UpdateCalfromUI, NICal_calstruct_init
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 18 July 2012 (SJS)
% 
% Revisions:
% 2 Aug 2012 (SJS): added CollectBackground to calibration settings
%--------------------------------------------------------------------------

%-----------------------------------------
% CALIBRATION SETTINGS
%-----------------------------------------
update_ui_val(handles.SideCtrl, cal.Side);
update_ui_str(handles.NrepsCtrl, cal.Nreps);
% Freq settings
update_ui_val(handles.FreqListCtrl, cal.UseFreqList);
if cal.UseFreqList
	disable_ui(handles.FminCtrl);
	disable_ui(handles.FmaxCtrl);
	disable_ui(handles.FstepCtrl);
else
	enable_ui(handles.FminCtrl);
	enable_ui(handles.FmaxCtrl);
	enable_ui(handles.FstepCtrl);
end
update_ui_str(handles.FminCtrl, cal.Fmin);
update_ui_str(handles.FmaxCtrl, cal.Fmax);
update_ui_str(handles.FstepCtrl, cal.Fstep);
% Attenuation settings
update_ui_val(handles.AttenFixCtrl, cal.AttenFix);
update_ui_str(handles.AttenFixValueCtrl, cal.AttenFixValue);
if cal.AttenFix
	show_uictrl(handles.AttenFixValueCtrl);
	set(handles.AttenFixValueCtrlText, 'Visible', 'on');
	disable_ui(handles.MinlevelCtrl);
	disable_ui(handles.MaxlevelCtrl);
	disable_ui(handles.AttenStepCtrl);
	disable_ui(handles.StartAttenCtrl);
else
	hide_uictrl(handles.AttenFixValueCtrl);
	set(handles.AttenFixValueCtrlText, 'Visible', 'off');
	enable_ui(handles.MinlevelCtrl);
	enable_ui(handles.MaxlevelCtrl);
	enable_ui(handles.AttenStepCtrl);
	enable_ui(handles.StartAttenCtrl);
end
update_ui_str(handles.MinlevelCtrl, cal.Minlevel);
update_ui_str(handles.MaxlevelCtrl, cal.Maxlevel);
update_ui_str(handles.AttenStepCtrl, cal.AttenStep);
update_ui_str(handles.StartAttenCtrl, cal.StartAtten);
%  Check Cal setting
update_ui_val(handles.CheckCalCtrl, cal.CheckCal + 1);
% Background
update_ui_val(handles.CollectBackgroundCtrl, cal.CollectBackground);
%-----------------------------------------
% INPUT/OUTPUT SETTINGS
%-----------------------------------------
update_ui_str(handles.ISICtrl, cal.ISI);
update_ui_str(handles.StimDurationCtrl, cal.StimDuration);
update_ui_str(handles.StimDelayCtrl, cal.StimDelay);
update_ui_str(handles.StimRampCtrl, cal.StimRamp);
update_ui_str(handles.SweepDurationCtrl, cal.SweepDuration);
update_ui_str(handles.DAscaleCtrl, cal.DAscale);
update_ui_val(handles.MeasureLeakCtrl, cal.MeasureLeak);
% input bandpass filter
update_ui_val(handles.InputFilterCtrl, cal.InputFilter);
update_ui_str(handles.HiPassFcCtrl, cal.InputHPFc);
update_ui_str(handles.LoPassFcCtrl, cal.InputLPFc);
if cal.InputFilter
	enable_ui(handles.HiPassFcCtrl);
	enable_ui(handles.LoPassFcCtrl);
else
	disable_ui(handles.HiPassFcCtrl);
	disable_ui(handles.LoPassFcCtrl);
end
%-----------------------------------------
% MICROPHONE SETTINGS
%-----------------------------------------
update_ui_val(handles.FRenableCtrl, cal.FRenable);
update_ui_val(handles.InputChannelCtrl, cal.InputChannel);
update_ui_str(handles.MicGainCtrl, cal.MicGain);
update_ui_str(handles.MicSensitivityCtrl, cal.MicSensitivity);		
%-----------------------------------------
% FILE SETTINGS
%-----------------------------------------
update_ui_str(handles.MicFRFileCtrl, cal.mic_fr_file);
update_ui_str(handles.CalFileCtrl, cal.calfile);
update_ui_val(handles.AutoSaveCtrl, cal.AutoSave);

