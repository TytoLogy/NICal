function Freqs = noctave(N, Fstart, Fend)
%------------------------------------------------------------------------
% Freqs = noctave
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% See also: DAQ Toolbox 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 3 August, 2012 (SJS)
%
% Revisions:
%------------------------------------------------------------------------


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% parse inputs
%------------------------------------------------------------------------
%------------------------------------------------------------------------
errStr = '';
if nargin ~= 3
	errStr = 'need N steps, Freq start and Freq end as inputs!';
end
if N < 1
	errStr = 'N must be greater than 0';
end
if Fstart > Fend
	errStr = 'Fstart must be less than Fend';
end
if Fstart <= 0
	errStr = 'Fstart must be greater than zero';
end
if Fend <= 0
	errStr = 'Fend must be greater than zero';
end

if ~isempty(errStr)
	error('%s: %s', mfilename, errStr);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% generate freqs
%------------------------------------------------------------------------
%------------------------------------------------------------------------

runFlag = 1;
f = 0;
n = 1;
while runFlag
	Fn = Fstart * power(2, f/N);
	if Fn <= Fend
		Freqs(n) = Fn;
		n = n + 1;
		f = f + 1;
	else
		runFlag = 0;
	end
	
	
end

