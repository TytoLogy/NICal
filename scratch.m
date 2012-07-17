

freqlist = [];

%--------------------------------------------------------------------------
% get filename and path to frequency list file (.m or .txt)
%--------------------------------------------------------------------------
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
	return
end

% load the numbers, store in tmp
try
	tmp = load(fullfile(pathname, filename));
catch errObj
	disp(errObj.message)
	return
end


freqlist = sort(unique(tmp))
nfreq = length(freqlist)


