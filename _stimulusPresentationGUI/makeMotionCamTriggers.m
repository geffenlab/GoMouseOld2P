% Make very large motion Camera trigger file
fs = 400000;
frameRate = 20;
pulse = [ones(0.001*fs,1)*3;zeros(0.049*fs,1)]; % 20 Hz frame rate
dur = round(75*60); % 90 mins in s
motionCamTriggers = repmat(pulse,frameRate*dur,1)/10;
audiowrite('temp_motionCamTrigger_400k.wav',motionCamTriggers,fs,...
    'bitsPerSample',16);