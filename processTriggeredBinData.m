function varargout = processTriggeredBinData(varargin)
%------------------------------------------------------------------------
% dbvals = processTriggeredBinData(	'filename', <file name for .bin file>,
%												'mode', <analysis mode>,
% 										'freqwidth' <automagic frequency detect width> )
%------------------------------------------------------------------------
% TytoLogy:NICal program
%------------------------------------------------------------------------
% FOR SESSION DATA ONLY
% 
% If inputfile is provided, program will read raw data from that file
% Otherwise, user will be prompted for input file
% 
% mode will select how the data are processed:
% 	mode = 'tones'  (default) will look for tone magnitudes in
% 					each .daq file and compute a freq-dB SPL curve
% 
% 	mode = 'rms'	will simply compute the overall db SPL level
% 							for each file
% 	mode = 'window' will compute db SPL levels for each file 
% 							broken up into 100 msec windows
% 
%	'rmswin', <value> where value is time window, in ms, to calculate 
%							running rms of signal
% 
% 	'stimchannel', <value>	default is 1
% 	'micchannel', <value>	default is 2
%
% freqwidth sets the width to search for peak frequency magnitude
% 		default is 21 Hz
%------------------------------------------------------------------------
% See also: NICal 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 14 April, 2017 from processTriggeredData (SJS)
%
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Variables and Constants Declarations
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% window size (in milliseconds) for computing rms (and dB) values
%	use smaller values for greater resolution, larger for coarse resolution
rms_windowsize_ms = 100;

% highpass cutoff frequency (Hz)
fcutoff = 90;
% filter order
forder = 3;

% Decimation factor - plotted data will be 1 / DeciFactor shorter
% and sampling rate will be Fs / DeciFactor
DeciFactor = 10;

% tolerance (in Hz) for finding peak frequency in autodetect mode
FreqDetectWidth = 21;

% stimulus channel
S = 1;
% mic channel
M = 2;

% set basepath and basename to empty and
% use default calibration mode ('tones')
inputfile = ''; 
calmode = 'tones';

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% parse inputs
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% loop through arguments
index = 1;
while index <= nargin
	% act according to tag
	switch lower(varargin{index})
		case 'inputfile'
			% set inputfile
			inputfile = varargin{index+1}; 
			% increment index by 2 places
			index = index + 2;
			
		case 'mode'
			switch(lower(varargin{index+1}))
				case {'tones', 'window', 'rms'}
					calmode = lower(varargin{index+1});
				otherwise
					error('%s: invalid mode value %s', mfilename, varargin{index+1});
			end
			index = index + 2;

		case 'freqwidth'
			FreqDetectWidth = varargin{index+1};
			index = index + 2;
			
		case 'rmswin'
			rms_windowsize_ms = varargin{index+1};
			index = index + 2;

		case 'stimchannel'
			S = varargin{index+1};
			index = index + 2;

		case 'micchannel'
			M = varargin{index+1};
			index = index + 2;

		otherwise
			error('%s: invalid option %s', mfilename, varargin{index});
	end	
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% find data file
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%----------------------------------
% check if filename was provided
%----------------------------------
if isempty(inputfile)
	% if empty, ask user for file
	%----------------------------------
	% output data path and file
	%----------------------------------
	DefaultPath = pwd;
	
	%----------------------------------
	% open panel to get .bin file name
	%----------------------------------
	[inputfile, basepath] = uigetfile(...
			 {'*.bin', 'BIN output files (*.bin)'}, ...
			  'Pick a .BIN file', DefaultPath);
	% check if user hit cancel (tmpfile, or basepath == 0)
	if isequal(inputfile, 0)
		disp('Cancelled...')
		varargout{1} = [];
		return
	else
		inputfile = fullfile(basepath, inputfile);
	end
	
else
	% inputfile was given, see if it exists 
	if ~exist(inputfile, 'file')
		error('%s: inputfile %s not found', mfilename, inputfile)
	end
end

%----------------------------------
% break down file into parts
%----------------------------------
[basepath, basename, ~] = fileparts(inputfile);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% check for matfile
%------------------------------------------------------------------------
%------------------------------------------------------------------------
matfile = fullfile(basepath, [basename '.mat']);
if ~exist(matfile, 'file')
	warning('%s: cannot locate NICal information file %s', mfilename, matfile);
	cal = [];
else
	load(matfile, '-MAT');
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% READ DATA
%------------------------------------------------------------------------
%------------------------------------------------------------------------
D = readBinData(inputfile);
[nSweeps, nChannels] = size(D.data);
% get sample rate from the cal struct
Fs = D.cal.Fs;
%------------------------------------------------------------------------
% get a highpass filter for processing the  data
%------------------------------------------------------------------------
% Nyquist frequency
fnyq = Fs/2;
% filter coefficients
[fcoeffb, fcoeffa] = butter(forder, fcutoff/fnyq, 'high');
% SPL conversion
VtoPa = 1 ./ (D.cal.Gain(1) * invdb(D.cal.MicGain(1)) * D.cal.MicSensitivity);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% should consider loading/asking user for 
% tone frequencies if tone is specified
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if strcmpi(calmode, 'tones')
	qVal = query_user('Auto-detect tone frequencies', 1);
	if qVal == 0
		% check if user has text list of frequencies
		fVal = query_user('Read list of frequencies from .txt file', 1);
		if fVal == 1
			% open panel to get .txt file name
			[freqfile, freqpath] = uigetfile(...
					 {'*.txt', 'txt frequency list files (*.txt)'}, ...
					  'Pick a .txt file');
			% check if user hit cancel (tmpfile, or basepath == 0)
			if isequal(freqfile, 0) || isequal(freqpath, 0)
				disp('Cancelled file load...')
				fVal = 0;
			else
				calfreqs = load(fullfile(freqpath, freqfile));
			end
		end
		% don't use else-if in order to allow fall-through from previous if
		if fVal == 0
			AUTOFREQ = 0;
			fstr = '';
			fstr = query_uservalue('Enter frequencies, separated by spaces', '');
			calfreqs = str2num(fstr); %#ok<ST2NM>
			clear fstr;			
		end
	else
		AUTOFREQ = 1;
		calfreqs = zeros(1, nSweeps);
	end
	if length(calfreqs) ~= nSweeps
		calfreqs = calfreqs * ones(1, nSweeps);
	end
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% PROCESS DATA
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%--------------------------------
% loop through daqfiles
%--------------------------------
for n = 1:nSweeps
	if nChannels == 2
		tmpdata = [D.data{n, 1} D.data{n, 2}];
	else
		tmpdata = [D.data{n, 1} D.data{n, 1}];
	end
	% ASSUME (!!) that stimulus data are collected on channel 1 (NI 0)
	% and mic data are on channel 2 (NI AI0)
	stimdata = tmpdata(:, S);
	micdata = tmpdata(:, M);
	clear tmpdata

	%------------------------------------------------------------------------
	% now process data
	%------------------------------------------------------------------------
	% window and filter the data
	micdata = sin2array(micdata', 1, Fs);
	micdata = filter(fcoeffb, fcoeffa, micdata);
	
	%--------------------------------
	% process data according to mode
	%--------------------------------
	switch lower(calmode)
		%--------------------------------
		% WINDOW
		%--------------------------------
		case 'window'
			[rms_vals, rmw_windows] =  processWindows(micdata, ...
																	rms_windowsize_ms, Fs);
			% convert to dB SPL
			dbvals{n} = dbspl(VtoPa * rms_vals); %#ok<AGROW>

			% decimate data for plotting
			micdata_reduced = decimate(micdata, DeciFactor);
			Fs_reduced = Fs / 10;
			% build time vectors for plotting
			t1 = ((1:length(micdata_reduced)) - 1) / Fs_reduced;
			t2 = rms_windowsize_ms * 0.001 * (0:size(dbvals{n},1)-1); %#ok<NASGU>
% 			t2 = rms_windowsize_ms * 0.001 * (0:rmsIndex);
% 			t2 = rms_windowsize_ms * 0.001 * (0:rmsIndex);
			t2 = rms_windowsize_ms * 0.001 * rmw_windows;
			% plot!
			figure(n)
			subplot(211)
			plot(t1, micdata_reduced);
			grid
			ylabel('Volts');
			subplot(212)
			plot(t2, dbvals{n}, 'Marker', '.', 'Color', 'r');
			ylabel('dB SPL')
			xlabel('Time (seconds)');
			grid
			title(sprintf('Peak dB SPL = %.2f', max(dbvals{n})))

		%--------------------------------
		% TONES
		%--------------------------------
		case 'tones'		
			[mags(n), phis(n), calfreqs(n)] = ...
				processTones(micdata, Fs, calfreqs(n), FreqDetectWidth); %#ok<AGROW>
			
		%--------------------------------
		% RMS
		%--------------------------------
		case 'rms'
			rms_vals(n) = rms(micdata);
			dbvals(n) = dbspl(VtoPa * rms_vals(n)); %#ok<AGROW>
		
	end
	
end


% assign output vars
switch lower(calmode)
	case 'tones'
		ndatums = length(calfreqs);
		out.freqs = calfreqs;
		out.mags = mags;
		out.phis = phis;
		out.dbvals = dbspl(VtoPa * rmssin * out.mags);

		figure
		subplot(211)
		plot(dbspl(VtoPa * rmssin * out.mags), '.-')
		set(gca, 'XTick', 1:ndatums);
		set(gca, 'XTickLabel', '');
		xlim([0 ndatums+1])
		grid
		title(basename);
		ylabel('dB SPL')

		subplot(212)
		plot(unwrap(out.phis), '.-')
		set(gca, 'XTick', 1:ndatums);
		labels = cell(ndatums, 1);
		for n = 1:ndatums
			labels{n} = sprintf('%.1f', 0.001*calfreqs(n));
		end
		set(gca, 'XTickLabel', labels);
		xlim([0 ndatums+1])
		grid
		ylabel('phase (rad)');
		xlabel('frequency (kHz)');

	case 'window'
		out.dbvals = dbvals;
		out.rms_windowsize_ms = rms_windowsize_ms;

	case 'rms'
		out.dbvals = dbvals;
		out.rms_vals = rms_vals;
end

out.Fs = Fs;
out.files = inputfile;
out.path = basepath;
out.fcutoff = fcutoff;
out.forder = forder;
out.DeciFactor = DeciFactor;
out.fcoeff.a = fcoeffa;
out.fcoeff.b = fcoeffb;
out.VtoPa = VtoPa;
out.data = D.data;
out.cal = D.cal;
varargout{1} = out;

end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% processTones function
%------------------------------------------------------------------------
%------------------------------------------------------------------------
function [mag, phi, freq] = processTones(micdata, Fs, calfreq, FreqDetectWidth)
	% get spectrum of data
	[tmpfreqs, tmpmags, fmax, magmax] = daqdbfft(micdata, Fs, length(micdata));
	
	if calfreq == 0
		% if  calfreq == 0, use fmax to detect magnitude (automatic freq. detection)
		freq = fmax;
	elseif FreqDetectWidth == 0
		freq = calfreq;
	else
		% otherwise, search in a range around the provided frequency
		freqindx = find(between(tmpfreqs, calfreq - FreqDetectWidth, calfreq + FreqDetectWidth));
		[~, maxindx] = max(tmpmags(freqindx));
		freq = tmpfreqs(freqindx(maxindx));
	end
	% compute pll magnitude and phase at this frequency
	[mag, phi] = fitsinvec(micdata, 1, Fs, freq);

end
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% processWindows function
%------------------------------------------------------------------------
%------------------------------------------------------------------------
function [rmsvals, rmswindows] = processWindows(data, windowms, fs)
	% convert windowsize from msec into samples
	windowpts = ms2samples(windowms, fs);
	
	% calculate rms_window indices into data
	rmswindows = 1:windowpts:length(data);
	Nwindows = length(rmswindows);
	
	% compute rms values for all windows of data
	rmsvals = zeros(Nwindows, 1);
	rmsIndex = 0;
	for w = 2:Nwindows
		rmsIndex = rmsIndex + 1;
		rmsvals(rmsIndex) = rms(data(rmswindows(w-1):rmswindows(w)));
	end
end
