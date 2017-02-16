function [start_bin, end_bin] = startendbins(handles)
%--------------------------------------------------------------------------
%  [start_bin, end_bin] = startendbins(handles)
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% computes start and end bins based on sample rate, handles.cal.StimDelay,
% handles.cal.StimRamp, handles.cal.StimDuration
%--------------------------------------------------------------------------
% See also: NICal
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	8 February, 2017
%
% Revisions:
%--------------------------------------------------------------------------

% start_bin depends on stimulus delay and onset ramp
start_bin = ms2bin(handles.cal.StimDelay + handles.cal.StimRamp, ...
								handles.iodev.Fs);
% if no start bin, set to 1
if ~start_bin
	start_bin = 1;
end
% end_bin depends on start bin and StimDuration (less Stimulus ramp)
end_bin = start_bin + ms2bin(handles.cal.StimDuration - handles.cal.StimRamp,...
											handles.iodev.Fs);
