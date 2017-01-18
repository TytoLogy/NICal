function [resp, index] = nidaq_calibration_io(iodev, stim, inpts)
%--------------------------------------------------------------------------
% [resp, index] = nidaq_calibration_io(iodev, stim, inpts)
%--------------------------------------------------------------------------
% TytoLogy -> Calibration -> NICal program
%--------------------------------------------------------------------------
% input/output for calibration
%------------------------------------------------------------------------
% Input Arguments:
%	iodev		input/output struct
%	stim		stimulus vector
%	inpts		# samples to collect (input)
% 
% Output Arguments:
%	resp		collected data {2X1} cell with vectors (1Xinpts) in size
%	index		# points collected
%------------------------------------------------------------------------
% See also: NICal, nidaq_aiao_init
%------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sharad J Shanbhag
% sshanbhag@neomed.edu
%--------------------------------------------------------------------------
% Created: 9 July 2012 (SJS)
% 
% Revisions:
%	18 Jan 2017 (SJS): updated comments
%--------------------------------------------------------------------------

% load stimulus onto NI memory
putdata(iodev.NI.ao, stim');
% calculate wait time from # of acquisition points
timeToWait = bin2seconds(inpts, iodev.Fs)*2;
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
