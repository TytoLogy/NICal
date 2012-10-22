function [DataPos, FreqList, calStruct] = readRawData(datfile)
%--------------------------------------------------------------------------
% readRawData.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 16 October, 2012
% 
% Revisions:
%	17 Oct 2012 (SJS): functionalized
%--------------------------------------------------------------------------

DataPos = [];
FreqList = [];
calStruct = [];

%-----------------------------------------------------------------------
% check dat file - if it exists, use it.  if not, check with user
%-----------------------------------------------------------------------
if ~exist(datfile, 'file')
	warning('NICAL:FileNotFound', '%s: could not find .DAT file', mfilename);
	fprintf('\t\tdat file: %s\n\n', datfile);

	[datfile, basepath] = uigetfile({'*.dat'},'Open dat file', ...
													basepath);
	
	if isequal(datfile, 0) || isequal(basepath, 0)
		return
	else
		[~, basename, ~] = fileparts(datfile);
		datfile = fullfile(basepath, [basename, '.dat']);
	end
else
	fprintf('Reading data information from %s\n', datfile);
end

%-----------------------------------------------------------------------
% read info from .dat file
%-----------------------------------------------------------------------
% open file
fp = fopen(datfile, 'r');
% read in cal struct (to get freqs and other info)
calStruct = readStruct(fp);
% build list of frequencies and file positions to data
% DataPos has dimensions [Nfreqs, Nreps]
FreqList = calStruct.Freqs;
DataPos = zeros(calStruct.Nfreqs, calStruct.Nreps);
for f = 1:calStruct.Nfreqs
	for r = 1:calStruct.Nreps
		% store position
		DataPos(f, r) = ftell(fp);
		% read data to move to next array
		tmp = readCell(fp);
	end
end
% close file
fclose(fp);
