function [freqlist, nfreqs, freqfile] = NICal_LoadFreqList(varargin)
%--------------------------------------------------------------------------
% NICal_LoadFreqList.m
%--------------------------------------------------------------------------
% NICal program
% TytoLogy Project
%--------------------------------------------------------------------------
% loads list of frequencies from text file
% if no filename is provided as input, will open a gui panel to ask
% user to select .m or .txt file that has a list of frequencies.
% 
% 	e.g.:
% 	
% 		20
% 		30
% 		40
% 		50
% 		60
% 		
%------------------------------------------------------------------------
% Input Arguments:
% 	freqfile		frequency text file (optional)
% 
% Output Arguments:
% 	freqlist		[nfreqs X 1] vector of unique frequencies from file,
% 					sorted from low to high
% 					
% 					[] if cancelled or error occurs
% 						
% 	nfreqs		# of unique frequencies from freqfile
% 	
% 	freqfile		filename (with path) of frequency list file
%------------------------------------------------------------------------
% See also: NICal
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 17 July, 2012 (SJS)
%
% Revisions:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% initialize vars
%--------------------------------------------------------------------------
freqlist = [];
nfreqs = 0;

%--------------------------------------------------------------------------
% check if user provided filename
%--------------------------------------------------------------------------
if ~nargin
	% get filename and path to frequency list file (.m or .txt)
	[filename, pathname] = uigetfile(	...
													{	...
														'*.txt', 'text files (*.txt)'; ...
														'*.m', 'Matlab m-files (*.m)'; ...
														'*.*', 'All Files (*.*)'	...
													}, ...
													'Pick a file', ...
													pwd	...
												);
	% return empty freqlist if user pressed "Cancel" button
	if any([(filename == 0), (pathname == 0)])
		freqfile = [];
		return
	else
		freqfile = fullfile(pathname, filename);
	end
	
else
	% use provided freqfile
	freqfile = varargin{1};
end

%--------------------------------------------------------------------------
% load the numbers, store in tmp (catch errors)
%--------------------------------------------------------------------------
try
	tmp = load(freqfile);
catch errObj
	disp(errObj.message)
	return
end

%--------------------------------------------------------------------------
% find unique values, sort from low to high
%--------------------------------------------------------------------------
freqlist = sort(unique(tmp));
nfreqs = length(freqlist);


