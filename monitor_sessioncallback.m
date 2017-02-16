function monitor_sessioncallback(obj, event) 
%------------------------------------------------------------------------
% monitor_sessioncallback(obj, event)
%------------------------------------------------------------------------
% TytoLogy -> NICal
%------------------------------------------------------------------------
% callback for monitoring continuous input
%------------------------------------------------------------------------
% See also: DAQ Toolbox 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 7 February, 2017 (SJS) from monitor_callback
%
% Revisions:
%------------------------------------------------------------------------

%---------------------------------------------------------------
% global variables
%---------------------------------------------------------------
global VtoPa Gain fcoeffa fcoeffb ...
		tvec_acq fvec Lacq Racq Lfft Rfft H SweepPoints %#ok<NUSED>

%---------------------------------------------------------------
% read data from ai object
%---------------------------------------------------------------
tmpdata = event.Data;
%---------------------------------------------------------------
% zero pad and filter data
%---------------------------------------------------------------
Lacq = buffer_filter(tmpdata(:, 1)', 5, obj.Rate, fcoeffb, fcoeffa);
Racq = buffer_filter(tmpdata(:, 2)', 5, obj.Rate, fcoeffb, fcoeffa);
%---------------------------------------------------------------
% find peak value
%---------------------------------------------------------------
maxval(1) = max(abs(Lacq));
maxval(2) = max(abs(Racq));
maxval = 1000 * maxval;
%---------------------------------------------------------------
% compute dB and update the dB SPL value
%---------------------------------------------------------------
dbSPLval(1) = dbspl(VtoPa * (rms(Lacq) ./ Gain(1)));
dbSPLval(2) = dbspl(VtoPa * (rms(Racq) ./ Gain(2)));
%---------------------------------------------------------------
% update data plot
%---------------------------------------------------------------
refreshdata(H.Lacq, 'caller');
refreshdata(H.Racq, 'caller');
%---------------------------------------------------------------
% plot fft
%---------------------------------------------------------------
[~, Lfft] = daqdbfft(Lacq, obj.Rate, SweepPoints);
[~, Rfft] = daqdbfft(Racq, obj.Rate, SweepPoints);
refreshdata(H.Lfft, 'caller');
refreshdata(H.Rfft, 'caller');
drawnow
%---------------------------------------------------------------
% update display on GUI
%---------------------------------------------------------------
update_ui_str(H.LValText, sprintf('%.4f', maxval(1)));
update_ui_str(H.LSPLText, sprintf('%.4f', dbSPLval(1)));
update_ui_str(H.RValText, sprintf('%.4f', maxval(2)));
update_ui_str(H.RSPLText, sprintf('%.4f', dbSPLval(2)));
