function testcallback(obj, event, lH) 
%------------------------------------------------------------------------
% testcallback(obj, event)
%------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal
%------------------------------------------------------------------------
% callback for logging triggered data
%
%	obj	caller (event source)
%	event	struct from caller (has Data)
%------------------------------------------------------------------------
% See also: NICal, NICal_RunTriggeredCalibration, DAQ Toolbox
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 7 February, 2017 (SJS) from monitor_sessioncallback
%
% Revisions:
%------------------------------------------------------------------------


%---------------------------------------------------------------
% read data from ai object and plot it
%---------------------------------------------------------------
fprintf('Triggered at %f\n', rem(event.TriggerTime, 1));

lH(1).XData = event.TimeStamps - event.TimeStamps(1);
lH(2).XData = lH(1).XData;
lH(1).YData = event.Data(:, 1);
lH(2).YData = event.Data(:, 2);

