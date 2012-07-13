function [resp, index] = nidaq_calibration_io(iodev, stim, inpts)
%--------------------------------------------------------------------------
% [resp, index] = nidaq_calibration_io(iodev, stim, inpts)
%--------------------------------------------------------------------------
% NICal program
% TytoLogy Project
% initializes nidaq system
%------------------------------------------------------------------------
% Input Arguments:
% 
% Output Arguments:
%------------------------------------------------------------------------
% See also: NICal
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 9 July 2012 (SJS)
% 
% Revisions:
%--------------------------------------------------------------------------

putdata(iodev.NI.ao, stim');

timeToWait = bin2seconds(inpts, iodev.Fs)*1.5;

%START ACQUIRING
start([iodev.NI.ai iodev.NI.ao]);
trigger([iodev.NI.ai iodev.NI.ao]);
wait(iodev.NI.ai, timeToWait);
% stop acquiring
stop([iodev.NI.ai iodev.NI.ao]);
% read data from ai object
index = iodev.NI.ai.SamplesAvailable;
rawdata = getdata(iodev.NI.ai, index);
resp{1} = rawdata(:, 1)';
resp{2} = rawdata(:, 2)';
