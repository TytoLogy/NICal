function plot_callback(obj, event)
%------------------------------------------------------------------------
% plot_callback
% 
%------------------------------------------------------------------------
% 
% 
%------------------------------------------------------------------------
% See also: DAQ Toolbox 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created:	5 July, 2012 (SJS)
%				modified from ai_plotpeek_callback.m
% Revisions:
%------------------------------------------------------------------------

%---------------------------------------------------------------
% global variables
%---------------------------------------------------------------
global tvec_acq Lacq Racq fvec Lfft Rfft H SweepPoints deciFactor start_bin end_bin

if strcmpi(obj.Running, 'On') || (obj.SamplesAvailable >= SweepPoints)
	%---------------------------------------------------------------
	% read data (acqpts, nchannels) array from ai object
	%---------------------------------------------------------------
	tmpdata = peekdata(obj, SweepPoints);

	%---------------------------------------------------------------
	% update data plot
	%---------------------------------------------------------------
	Lacq = downsample(tmpdata(:, 1), deciFactor);
	Racq = downsample(tmpdata(:, 2), deciFactor);

	refreshdata(H.Lacq, 'caller');
	refreshdata(H.Racq, 'caller');
	
	% plot fft
	[tmpf, Lfft] = daqdbfft(tmpdata(start_bin:end_bin, 1), obj.SampleRate, length(start_bin:end_bin));
	[tmpf, Rfft] = daqdbfft(tmpdata(start_bin:end_bin, 2), obj.SampleRate, length(start_bin:end_bin));
	refreshdata(H.Lfft, 'caller');
	refreshdata(H.Rfft, 'caller');
	drawnow

end