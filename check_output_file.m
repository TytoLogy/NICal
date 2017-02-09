function calfile = check_output_file(handles)
%--------------------------------------------------------------------------
% calfile = check_output_file(handles)
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% checks calibration file in handles.cal.calfile
%	if it exists, asks if necessary to overwrite
%			if new file, get one from user
% returns calfile OR 0 if cancelled
%--------------------------------------------------------------------------
% See also: NICal
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created:	8 February, 2017
%
% Revisions:
%--------------------------------------------------------------------------

% check if file exists
if exist(handles.cal.calfile, 'file')
	% if it exists, ask if necessary to overwrite
	resp = uiyesno('title', 'Save File', 'string', ...
							'File exists! Overwrite?', 'default', 'No');
	% if no, get a new file
	if strcmpi(resp, 'No')
		[pathstr, fname, fext] = fileparts(handles.cal.calfile);
		[newname, newpath] = uiputfile('*.cal', ...
													'Save calibration data to file', ...
													fullfile(pathstr, [fname '_1' fext]));
		% check for cancellation
		if isequal(newname, 0) || isequal(newpath, 0)
			calfile = 0;
			return
		else
			calfile = fullfile(newpath, newname);
		end
	end
else
	calfile = handles.cal.calfile;
end
