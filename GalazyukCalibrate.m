
% 15 Nov 2017 (SJS): updating to use processTriggeredBinData 
% caldata=processTriggeredData;
caldata=processTriggeredBinData;

for i=1:length(caldata.freqs)
	Freq_list(i)=caldata.freqs(i)/97656.25;
end

Freq_list(1)=0;
Freq_list(93)=1;
Freq_list=Freq_list';


Gain_list=caldata.mags;
Gain_list(1)=Gain_list(2);
Gain_list(93)=Gain_list(92);
Gain_list=Gain_list';

ntaps= 250;

nyquist= 97656.25;

filtcoefs = fir2(ntaps, Freq_list, 10.^(Gain_list/20));

fileid = fopen('D:\Galazyuk\GalazyukCalibration.txt','wt+');

fprintf(fileid, '%6f\n', filtcoefs);

fclose(fileid);

