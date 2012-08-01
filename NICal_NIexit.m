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
EventLogAI = showdaqevents(handles.iodev.NI.ai);
EventLogAO = showdaqevents(handles.iodev.NI.ao);

% delete and clear ai and ch0 object
delete(handles.iodev.NI.ai);
delete(handles.iodev.NI.ao);
delete(handles.iodev.NI.chI);
delete(handles.iodev.NI.chO);
clear handles.iodev.NI.ai handles.iodev.NI.ao handles.iodev.NI.chI handles.iodev.NI.chO

% save settings information to mat file
save(fullfile(pwd, 'NICal_EventLogs.mat'), ...
		'EventLogAI'			, ...
		'EventLogAO'			, ...
		'-MAT' );


	
