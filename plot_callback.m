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
global Lacq Racq H acqpts

if strcmpi(obj.Running, 'On')
	%---------------------------------------------------------------
	% read data from ai object
	%---------------------------------------------------------------
	tmpdata = peekdata(obj, acqpts);
	data{1} = tmpdata(:, 1);
	data{2} = tmpdata(:, 2);
	
	%---------------------------------------------------------------
	% update data plot
	%---------------------------------------------------------------
	refreshdata(Lacq);
	refreshdata(Racq);
	
end