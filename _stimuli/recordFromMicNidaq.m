%% setup nidaq

clear all

s = daq.createSession('ni');
fs = 400000;
s.Rate = fs;
ao1=addAnalogOutputChannel(s,'dev1',0,'Voltage');
% ao2=addAnalogOutputChannel(s,'dev1',1,'Voltage');
ai1=addAnalogInputChannel(s,'dev1',1,'Voltage');
% ai2=addAnalogInputChannel(s,'dev1',8,'Voltage');

s.Channels(2).InputType = 'SingleEnded';
% s.Channels(4).InputType = 'SingleEnded';


ref_Pa = 20e-6;
volts_per_Pa = .316;
%% Record
% stim=audioread('C:\Users\Chris\Documents\GitHub\Kath\stimuli\randReg_pool580dB_1rep.wav');
% stim=stim*10;
output = [zeros(fs*10,1)];
% output = stim(1:fs*10,1);
queueOutputData(s,[output]);
[data,startTime] = s.startForeground();
[fb, fa] = butter(5, 2*300 / fs, 'high');
b = filter(fb, fa, data) / ref_Pa / volts_per_Pa;
b = b - mean(b);
% noise_ms = mean(b.^2);


%% analyse

% repsM = data/ref_PA/volts_per_PA;
[P,f] = pwelch(b, 1024, 120, [], fs, 'onesided');
dB = 10*log10(P);
plot(f,dB,'r');
title(['Total volume ' num2str(10*log10(mean(P)*(f(end)-f(1))))...
    'dB in response to background noise with galvo-galvo scanning.']);




