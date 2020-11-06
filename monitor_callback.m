function monitor_callback(obj, event)
%------------------------------------------------------------------------
%  monitor_callback(obj, event)
%------------------------------------------------------------------------
% TytoLogy -> NICal
%------------------------------------------------------------------------
% callback for monitoring continuous input% 
%------------------------------------------------------------------------
% See also: DAQ Toolbox 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 28 August, 2012 (SJS) from ai_plot_callback
%
% Revisions:
%	13 Sep 2012 (SJS): fixed incorrect Rfft calculation from daqdbfft
%	23 Mar 2017 (SJS): fixed filter on/off issue
%------------------------------------------------------------------------

% global variables
global VtoPa Gain fcoeffa fcoeffb filtEnable ...
		tvec_acq fvec Lacq Racq Lfft Rfft H SweepPoints

% read data from ai object
tmpdata = getdata(obj, SweepPoints);
% Lacq = filtfilt(fcoeffb, fcoeffa, sin2array(tmpdata(:, 1)', 5, obj.SampleRate));
% Racq = filtfilt(fcoeffb, fcoeffa, sin2array(tmpdata(:, 2)', 5, obj.SampleRate));

if filtEnable
	Lacq = buffer_filter(tmpdata(:, 1)', 5, obj.SampleRate, fcoeffb, fcoeffa);
	Racq = buffer_filter(tmpdata(:, 2)', 5, obj.SampleRate, fcoeffb, fcoeffa);
else
	Lacq = tmpdata(:, 1)';
	Racq = tmpdata(:, 2)';
end

% find peak value
maxval(1) = max(abs(Lacq));
maxval(2) = max(abs(Racq));
maxval = 1000 * maxval;
% compute dB and update the dB SPL value
dbSPLval(1) = dbspl(VtoPa * (rms(Lacq) ./ Gain(1)));
dbSPLval(2) = dbspl(VtoPa * (rms(Racq) ./ Gain(2)));

%---------------------------------------------------------------
% update data plot
%---------------------------------------------------------------
refreshdata(H.Lacq, 'caller');
refreshdata(H.Racq, 'caller');
% plot fft
[tmpf, Lfft] = daqdbfft(Lacq, obj.SampleRate, SweepPoints);
[tmpf, Rfft] = daqdbfft(Racq, obj.SampleRate, SweepPoints);
refreshdata(H.Lfft, 'caller');
refreshdata(H.Rfft, 'caller');
drawnow

update_ui_str(H.LValText, sprintf('%.4f', maxval(1)));
update_ui_str(H.LSPLText, sprintf('%.4f', dbSPLval(1)));
update_ui_str(H.RValText, sprintf('%.4f', maxval(2)));
update_ui_str(H.RSPLText, sprintf('%.4f', dbSPLval(2)));
