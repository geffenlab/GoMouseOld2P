function output = makeNoiseClick(dur,ISI,fs,FILT)

clickDur = 0.02; % 20 ms
clickRate = 10; % Hz

noise = randn(1,fs)/sqrt(2);
noise = (noise - mean(noise));
noise = noise.*10^(-(-10/20));
noise = noise(1,0.5*fs:(0.5+clickDur)*fs);
noise = envelopeKCW(noise(1:0.02*fs),5,fs);
noise = conv(noise,FILT,'same');
noise = [noise zeros(1,((1/clickRate)-clickDur)*fs)];
ch1 = [repmat(noise,1,clickRate*dur) zeros(1,fs*ISI)];
ch2 = [ones(1,dur*fs)*5,zeros(1,ISI*fs)];

output = [ch1;ch2]';
