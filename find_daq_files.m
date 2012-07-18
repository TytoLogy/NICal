function varargout = find_daq_files(basename, varargin)
%------------------------------------------------------------------------
% [daqFileNames, daqFileNumbers] = find_daq_files(basename, basepath)
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
% Created: 6 July, 2012 (SJS)
%
% Revisions:
%------------------------------------------------------------------------

if nargin == 2
	basepath = varargin{1};
else
	basepath = pwd;
end

% build search string for dir command from the basename
searchStr = sprintf('%s-*.daq', basename);
% get struct array from dir command
tmpstruct = dir(fullfile(basepath, searchStr));
% convert to cell array
tmpcell = struct2cell(tmpstruct);
% keep the names
daqFileNames = tmpcell(1, :);
% clear unused variables
clear tmpstruct tmpcell;

% get # of found daqFileNames
nDaqFiles = length(daqFileNames);
if nDaqFiles == 0
	varargout{1} = '';
	varargout{2} = [];
	return
end

% allocate vector for daq file numbers
daqFileNumbers = zeros(nDaqFiles, 1);
% loop through the daq file names
for n = 1:nDaqFiles
	% get the filename without the extension
	[~, tmpname, ~] = fileparts(daqFileNames{n});
	% find the dashes in the filename
	dashindices = find(tmpname == '-');
	% take the characters after the final dash and convert to number,
	% then store this number in daqFileNumbers
	daqFileNumbers(n) = str2num(tmpname( (dashindices(end)+1):end));
end

% sort the numbers and reorder the daqFileNames
[sorted_vals, sorted_indices] = sort(daqFileNumbers);
varargout{2} = daqFileNumbers(sorted_indices);
varargout{1} = daqFileNames(sorted_indices);







