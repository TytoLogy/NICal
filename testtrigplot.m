%c = buildCalFromTriggeredData('NICaldata-19-Sep-2012', '/Users/sshanbhag/Work/Data/Audio/BatRig/19Sep2012/TDT3981_5Vpk_AvisoftMic3Chan0/');

xrange = 0.001*[min(c.freqs) max(c.freqs)];
VtoPa = 1/c.calibration_settings.MicSensitivity;

figure(2)
subplot(211)
plot(0.001*c.freqs, dbspl(VtoPa * rmssin * c.mags), '.-')
set(gca, 'XTickLabel', '');
grid
ylabel('dB SPL')
xlim(xrange);
legend('Ch0', 'Ch1')


subplot(212)
plot(0.001*c.freqs, unwrap(c.phis), '.-')
xlim(xrange);
grid
ylabel('phase (rad)');
xlabel('frequency (kHz)');


figure(3)
subplot(211)
plot(0.001*c.freqs, dbspl(VtoPa * rmssin * c.dist), '.-')
set(gca, 'XTickLabel', '');
grid
ylabel('dB SPL')
xlim(xrange);
legend('Ch0', 'Ch1')


subplot(212)
plot(0.001*c.freqs,  dbspl(VtoPa * rmssin * c.dist3), '.-')
xlim(xrange);
grid
xlabel('frequency (kHz)');
