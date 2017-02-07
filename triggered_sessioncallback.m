function triggered_sessioncallback(obj, event, datafile) 
%------------------------------------------------------------------------
% triggered_sessioncallback(obj, event, fid)
%------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal
%------------------------------------------------------------------------
% callback for logging triggered data
%
%	obj	caller (event source)
%	event	struct from caller (has Data)
%	T		struct with fields:
%				VtoPa		conversion factor from volts to Pascal
% 				Gain			input gain
% 				fcoeffa		filter coefficients
% 				fcoeffb		filter coefficients
% 				tvec_acq	time vector for plots
% 				fvec			frequency vector for plots
% 				Lacq			L channel acquired data in plot
% 				Racq			R channel acquired data in plot
% 				Lfft			L channel FFT data in plot
% 				Rfft			R channel FFT data in plot
% 				H				struct with handles Lacq and Racq for plots
%				side			channel being calibrated
%									1 = L, 2 = R, 3 = Both
% 				SweepPoints # of points per "sweep"
%	datafile	logging data file
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
% global variables
%---------------------------------------------------------------
 global VtoPa Gain fcoeffa fcoeffb ...
			Lacq Racq Lfft Rfft H SweepPoints side

%---------------------------------------------------------------
% read data from ai object and write to file
%---------------------------------------------------------------
tmpdata = event.Data;
% open file
fp = fopen(datafile, 'a');
% write L data
if any(side == [1 3])
	writeVector(fp, tmpdata(:, 1), 'double');
end
% write R data
if any(side == [2 3])
	writeVector(fp, tmpdata(:, 2), 'double');
end
% close file
fclose(fp);

%---------------------------------------------------------------
% Update plots with new data
%---------------------------------------------------------------
if any(side == [1 3])
	% zero pad and filter data
	Lacq = buffer_filter(tmpdata(:, 1)', 5, obj.Rate, fcoeffb, fcoeffa);
	% find peak value (mV)
	maxval(1) = 1000*max(abs(Lacq));
	% compute dB and update the dB SPL value
	dbSPLval(1) = dbspl(VtoPa * (rms(Lacq) ./ Gain(1)));
	% update fft data
	[~, Lfft] = daqdbfft(Lacq, obj.Rate, SweepPoints);
	% update display on GUI
	update_ui_str(H.LValText, sprintf('%.4f', maxval(1)));
	update_ui_str(H.LSPLText, sprintf('%.2f', dbSPLval(1)));
end
if any(side == [2 3])
	% zero pad and filter data
	Racq = buffer_filter(tmpdata(:, 2)', 5, obj.Rate, fcoeffb, fcoeffa);
	% find peak value (mV)
	maxval(2) = 1000*max(abs(Racq));
	% compute dB and update the dB SPL value
	dbSPLval(2) = dbspl(VtoPa * (rms(Racq) ./ Gain(2)));
	% update fft data
	[~, Rfft] = daqdbfft(Racq, obj.Rate, SweepPoints);
	% update display on GUI
	update_ui_str(H.RValText, sprintf('%.4f', maxval(2)));
	update_ui_str(H.RSPLText, sprintf('%.2f', dbSPLval(2)));
end
%---------------------------------------------------------------
% refresh data plot
%---------------------------------------------------------------
refreshdata(H.Lacq, 'caller');
refreshdata(H.Racq, 'caller');
%---------------------------------------------------------------
% refresh fft plot
%---------------------------------------------------------------
refreshdata(H.Lfft, 'caller');
refreshdata(H.Rfft, 'caller');
drawnow

