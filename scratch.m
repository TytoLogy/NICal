cal = cal_structinit;


s = normalize(synmononoise_fft(500, cal.Fs, 1, cal.Fs/2, 1, 0));
s = sin2array(s, 1, cal.Fs);

% filter coefficients using butterworth bandpass filter
S = filtfilt(cal.fcoeffb, cal.fcoeffa, s);

t = ((1:length(s)) - 1) / cal.Fs;

figure(1)
plot(t, s, t, S)
legend({'s', 'S'});
fftplot(s, cal.Fs, figure(2));
fftplot(S, cal.Fs, figure(3));
