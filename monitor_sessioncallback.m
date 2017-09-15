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
%	23 Mar 2017 (SJS): fixed filter on/off issue
%------------------------------------------------------------------------

%---------------------------------------------------------------
% global variables
%---------------------------------------------------------------
global VtoPa Gain fcoeffa fcoeffb filtEnable ...
		tvec_acq fvec Lacq Racq Lfft Rfft H SweepPoints ...
		SaveSweep SweepFile %#ok<NUSED>

%---------------------------------------------------------------
% read data from ai object
%---------------------------------------------------------------
tmpdata = event.Data;
%---------------------------------------------------------------
% zero pad and filter data
%---------------------------------------------------------------
if filtEnable
	Lacq = buffer_filter(sin2array(tmpdata(:, 1)', 2, obj.Rate), ...
														5, obj.Rate, fcoeffb, fcoeffa);
	Racq = buffer_filter(sin2array(tmpdata(:, 2)', 2, obj.Rate), ...
														5, obj.Rate, fcoeffb, fcoeffa);
else
	Lacq = tmpdata(:, 1)';
	Racq = tmpdata(:, 2)';
end
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
%---------------------------------------------------------------
% save sweep?
%---------------------------------------------------------------
if SaveSweep
	save(SweepFile, 'Lacq', 'Racq', 'tvec_acq', 'Lfft', 'Rfft', 'fvec', '-MAT');
end

