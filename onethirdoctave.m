function Freqs = onethirdoctave(Fstart, Fend)
%------------------------------------------------------------------------
% Freqs = onethirdoctave
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% See also: noctave
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
if nargin ~= 2
	errStr = 'need Freq start and Freq end as inputs!';
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

Freqs = noctave(3, Fstart, Fend);
% 
% runFlag = 1;
% f = 0;
% n = 1;
% while runFlag
% 	Fn = Fstart * power(2, f/3);
% 	if Fn <= Fend
% 		Freqs(n) = Fn;
% 		n = n + 1;
% 		f = f + 1;
% 	else
% 		runFlag = 0;
% 	end
% 	
% 	
% end
% 
% varargout{1} = Freqs;
