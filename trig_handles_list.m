%------------------------------------------------------------------------
% trig_handles_list
%------------------------------------------------------------------------
% NICal
% 
%------------------------------------------------------------------------
% See also: DAQ Toolbox 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: ??? 2012 (SJS)
%
% Revisions:
%	29 Nov 2012 (SJS): added some new fields, added comments
%------------------------------------------------------------------------
trigger_handle_names = {	...
	'AttenFixCtrl', ...
	'AttenFixValueCtrl', ...
	'AttenFixValueCtrlText', ...
	'AttenStepCtrl', ...
	'StartAttenCtrl', ...
	'DAscaleCtrl', ...
	'FmaxCtrl', ...
	'FminCtrl', ...
	'FstepCtrl', ...
	'FreqListCtrl', ...
	'FreqValText', ...
	'FstepCtrl', ...
	'ISICtrl', ...
	'MinlevelCtrl', ...
	'MaxlevelCtrl', ...	
	'NrepsCtrl', ...
	'SideCtrl', ...
	'CollectBackgroundCtrl', ...
	'MeasureLeakCtrl', ...
	'SaveRawDataCtrl', ...
};


trigger_handles = [	...
	handles.AttenFixCtrl, ...
	handles.AttenFixValueCtrl, ...
	handles.AttenFixValueCtrlText, ...
	handles.AttenStepCtrl, ...
	handles.StartAttenCtrl, ...
	handles.DAscaleCtrl, ...
	handles.FmaxCtrl, ...
	handles.FminCtrl, ...
	handles.FstepCtrl, ...
	handles.FreqListCtrl, ...
	handles.FreqValText, ...
	handles.FstepCtrl, ...
	handles.ISICtrl, ...
	handles.MinlevelCtrl, ...
	handles.MaxlevelCtrl, ...
	handles.NrepsCtrl, ...
	handles.SideCtrl, ...
	handles.CollectBackgroundCtrl, ...
	handles.MeasureLeakCtrl, ...
	handles.SaveRawDataCtrl, ...
];