% ask user for input file name
% default path
binpath = 'D:\Galazyuk';
binfile = [];
[binfile, binpath] = uigetfile('*.bin', ...
											'Read triggered data from file', ...
											fullfile(binpath, binfile));
if binfile == 0
	fprintf('read bin data cancelled\n');
	return
else
	bindatafile = fullfile(binpath, binfile);
end


% 15 Nov 2017 (SJS): updating to use processTriggeredBinData 
% caldata=processTriggeredData;
caldata=processTriggeredBinData('inputfile', bindatafile);

for i=1:length(caldata.freqs)
	Freq_list(i)=caldata.freqs(i)/97656.25;
end

% Assign frequencies and gain values
Freq_list(1)=0;
Freq_list(93)=1;
Freq_list=force_row(Freq_list);
Gain_list=caldata.dbvals;
Gain_list(1)=Gain_list(2);
Gain_list(93)=Gain_list(92);
Gain_list=force_row(Gain_list);
% plot
figure
subplot(211)
plot(Freq_list, Gain_list, '.-');
grid('on');
xlabel('Frequency (normalized)');
ylabel('dB SPL');

% convert gain to correction (attenuation values)
% first, find minimum value
mindbspl = min(Gain_list);
% correction = mindbspl - Gain
Gain_correction = mindbspl - Gain_list;
subplot(212)
plot(Freq_list, Gain_correction, '.-');
grid('on');
xlabel('Frequency (normalized)');
ylabel('correction dB');


% now build FIR filter to compensate
ntaps= 250;
% this is from TDT RZ6 hardware
nyquist= 97656.25;

filtcoefs = fir2(ntaps, Freq_list, 10.^(Gain_correction/20));

% ask user for output file name
% default path
outpath = 'D:\Galazyuk';
outfile = 'GalazyukCalibration.txt';
[outfile, outpath] = uiputfile('*.txt', ...
											'Save Filter Coefficients to File', ...
											fullfile(outpath, outfile));
if outfile == 0
	fprintf('Save filter coefficients cancelled\n');
	return
end
% write to coefficient file
fileid = fopen(fullfile(outpath, outfile),'wt+');
fprintf(fileid, '%6f\n', filtcoefs);
fclose(fileid);

