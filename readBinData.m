function D = readBinData(inputfile)
%------------------------------------------------------------------------
% D = readBinData(<input .bin file>)
%------------------------------------------------------------------------
% TytoLogy:NICal program
%------------------------------------------------------------------------
% FOR SESSION DATA ONLY!!!!!!!
% 
% Program will read raw data from that file
% 
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

if isempty(inputfile) || ~exist(inputfile, 'file')
	error('%s: empty or nonexisting file %s', mfilename, inputfile)
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% READ DATA
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%-----------------------------------------------
% open file
%-----------------------------------------------
fp = fopen(inputfile, 'r');

%-----------------------------------------------
% read cal struct
%-----------------------------------------------
tmp = readStruct(fp);
cal = tmp;
% convert character fields
cal.mic_fr_file = char(cal.mic_fr_file);
cal.calfile = char(cal.calfile);
cal.TriggerSettings.TriggerType = char(cal.TriggerSettings.TriggerType);
cal.TriggerSettings.TriggerSource = char(cal.TriggerSettings.TriggerSource);
cal.TriggerSettings.TriggerCondition = ...
									char(cal.TriggerSettings.TriggerCondition);
clear tmp

%-----------------------------------------------
% read sweeps
%-----------------------------------------------
% allocate temporary storage
if any(cal.Side == [1 2])
	tmp = cell(1000, 1);
else
	tmp = cell(1000, 2);	
end

fFlag = 0;
index = 1;
while(~feof(fp) && ~fFlag)
	% read L data
	if feof(fp)
		fFlag = 1;
	else
		if any(cal.Side == [1 3])
			tmp{index, 1} = readVector(fp);
		end
	end
	
	% read R data
	if feof(fp)
		fFlag = 1;
	else
		if any(cal.Side == [2 3])
			tmp{index, 2} = readVector(fp);
		end
	end
	
	% check for eof or error
	[fmsg, ferr] = ferror(fp);
	if feof(fp) || ferr
		fFlag = 1;
		fprintf('%s: end of file... message: <%s>\n', mfilename, fmsg);
	else
		index = index + 1;
	end
end

%-----------------------------------------------
% close file
%-----------------------------------------------
fclose(fp);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% assign to D struct
%------------------------------------------------------------------------
%------------------------------------------------------------------------
D.cal = cal;
% check that last sweep is not empty
% allocate temporary storage
if cal.Side == 3
	if (isempty(tmp{index, 1}) && isempty(tmp{index, 2}))
		% if so, move back 1
		index = index - 1;
		if index < 1
			error('%s: odd, or empty data file %s', mfilename, inputfile);
		end
	end
elseif cal.Side == 1
	if isempty(tmp{index, 1})
		% if so, move back 1
		index = index - 1;
		if index < 1
			error('%s: odd, or empty data file %s', mfilename, inputfile);
		end
	end
elseif cal.Side == 2
	if isempty(tmp{index, 2})
		% if so, move back 1
		index = index - 1;
		if index < 1
			error('%s: odd, or empty data file %s', mfilename, inputfile);
		end
	end
end
% store only valid data
D.data = tmp(1:index, :);
