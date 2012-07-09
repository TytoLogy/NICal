%--------------------------------------------------------------------------
% NICal_NIexit.m
%--------------------------------------------------------------------------
%
% closes NI devices nicely
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 9 July 2012 (SJS)
% 				Created from SpeakerCal_tdtexit.m
% 
% Revisions:
%	9 July, 2012 (SJS) renamed for NICal project
%--------------------------------------------------------------------------



%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% Clean up the RP circuits
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
disp('...closing NI devices...');


% get event log
EventLogAI = showdaqevents(iodev.NI.ai);
EventLogAO = showdaqevents(iodev.NI.ao);

% delete and clear ai and ch0 object
delete(iodev.NI.ai);
delete(iodev.NI.ao);
delete(iodev.NI.chI);
delete(iodev.NI.chO);
clear iodev.NI.ai iodev.NI.ao iodev.NI.chI iodev.NI.chO

% save settings information to mat file
save(fullfile(pwd, 'NICal_EventLogs.mat'), ...
		'EventLogAI'			, ...
		'EventLogAO'			, ...
		'-MAT' );


	
