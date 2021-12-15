% testMouseVocs
close all
clear all

fs = 400e3;
s = startSession(fs);
fs = s.Rate;
load('C:\calibration\Filters\20170420_2PspkrNidaqInvFilt_3k-80k_fs400k.mat');
n = 10;
targetVol = 70;
upperFreq = 70e3;
lowerFreq = 3000;
softGain = 10;
ref_PA = 20e-6;
volts_per_PA = .316;

% load mouse voc
fileLoc = 'C:\Users\Chris\Dropbox\GeffenLab\Kath\mouseVocalisations\extracted\mats\';
vfiles = dir([fileLoc '*.mat']);
load([fileLoc vfiles(5).name]);
voc = (voc-mean(voc))/max(abs(voc));
figure; spectrogram(voc,200,20,1000,fs,'yaxis');
title(['purified vocalisation'])
savefig('exampleMouseVoc005_purified.fig')

% RECORD SOME BACKGROUND NOISE
[fb, fa] = butter(5, 2*300 / fs, 'high');
disp('Acquiring 3s of background noise:');
silence = zeros(3*fs,1);
queueOutputData(s,silence);
[S, time] = s.startForeground();
b = filter(fb, fa, S) / ref_PA / volts_per_PA; % convert from Voltage to pressure
b = b - mean(b);
noise_ms = mean(b.^2);

% Record the mouse voc w/o filter
queueOutputData(s,voc);
[S, time] = s.startForeground();
[fb, fa] = butter(5, 2*300 / fs, 'high');
S = filter(fb,fa,S);
% t = S(4000:30000);
S = S-mean(S);
tt = S/ref_PA/volts_per_PA;
r = sqrt(mean((tt).^2)- noise_ms);
dB = 20*log10(r);
figure; spectrogram(S,200,20,1000,fs,'yaxis');
title(['Unfiltered: ' num2str(dB) ' dB'])
% repsM = S/ref_PA/volts_per_PA;
% [P,f] = pwelch(repsM, 1024, 120, [], fs, 'onesided');
% dB = 10*log10(mean(P(P>70)));
savefig('exampleMouseVoc005_unfiltered.fig')

% Record mouse voc filtered
vocf = conv(voc,FILT,'same');
vocf = vocf/0.5;
queueOutputData(s,vocf);
[S, time] = s.startForeground();
[fb, fa] = butter(5, 2*300 / fs, 'high');
S = filter(fb,fa,S);
S = S-mean(S);
tt = S/ref_PA/volts_per_PA;
% r = rms(tt);
r = sqrt(mean((tt).^2)- noise_ms);
dB = 20*log10(r)
figure; spectrogram(S,200,20,1000,fs,'yaxis');
title(['Filtered: ' num2str(dB) ' dB'])
savefig('exampleMouseVoc005_filtered.fig')

% Record mouse voc filtered attenuated
attn = -15;
voc = voc.*10^(-attn/20); % attenuate
vocf = conv(voc,FILT,'same');
vocf = vocf/0.5;
queueOutputData(s,vocf);
[S, time] = s.startForeground();
[fb, fa] = butter(5, 2*300 / fs, 'high');
S = filter(fb,fa,S);
S = S-mean(S);
tt = S/ref_PA/volts_per_PA;
% r = rms(tt);
r = sqrt(mean((tt).^2)- noise_ms);
dB = 20*log10(r)
figure; spectrogram(S,200,20,1000,fs,'yaxis');
title(['Filtered: ' num2str(dB) ' dB'])
savefig('exampleMouseVoc005_filtered&attenuated.fig')




