clear all

fs = 500000;
NFFT = 1024;
T = [18 1];



t = 0:1/fs:.250;
s = chirp(t, 1, max(t), fs/2);
N  = length(s); 
 
figure(1)
subplot(211)
tic
a = myspectrogram(s, fs, T, @hamming, NFFT, [-100 -1], false, 'default', false, 'per');
toc

subplot(212)
Tw = T(1);                          % frame width (ms)
Ts = T(2);                          % frame shift (ms)
Nw = round(fs*Tw*0.001);            % frame width (samples)
Ns = round(fs*Ts*0.001);            % frame shift (samples)
w = hamming(Nw).';
[S,F,T] = spectrogram(s, w(Nw).',Nw-Ns, NFFT,fs);
S = abs(S);                         % compute magnitude spectrum 
S = S/max(max(S));                  % normalize magntide spectrum
S = 20*log10(S);                    % compute power spectrum in dB


handle = imagesc(T, F, S, [Smin Smax]);
%handle = imagesc(T, F, S, 'CDataMapping', 'direct');
axis('xy');
axis([0 N/fs  0 fs/2]);
