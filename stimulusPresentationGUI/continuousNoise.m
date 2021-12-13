function continuousNoise(src,event,noiseClickDur,ISI,fs,FILT)

noise = makeNoiseClick(noiseClickDur,ISI,fs,FILT); % duration, ISI and sample rate
src.queueOutputData(noise);