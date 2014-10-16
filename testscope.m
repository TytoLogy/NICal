clear all

fs = 500000;
NFFT = 1024;
T = [10 5];



t = 0:1/fs:.250;
s = chirp(t, 1, max(t), fs/2);
N  = length(s); 
 
figure(1)
tic
a = myspectrogram(s, fs, T, @hamming, NFFT, [-100 -1], false, 'default', false, 'per');
toc

