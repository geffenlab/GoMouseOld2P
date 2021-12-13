% Generate stimuli for rate level function

clear
seed = 5;
rng(seed);

filename = 'noiseClicks_1rep_10clicks7sisi100dB_400k';
fs = 400000;
nDur = 25; % noise duration in ms
nClicks = 10; % number of noise clicks
clickRate = 10; % in Hz
ISI = 3000; % inter stim interval duration in ms
repeats = 1; % number of repeats of each stim
intensity = 100; % dB

% Load the calibration filter
filterName = 'C:\calibration\Filters\20170216_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(filterName);

% Create stimuli
totalDur = ((1/clickRate)*nClicks)+ISI/1000; % total duration in s

% preallocate variable
stimArray = zeros(1,round(totalDur*fs));
events = stimArray;
loc = 1; % placement of the tone in the zeros
n = rand(1,(nDur/1000)*fs);
n = envelopeKCW(n,5,fs); % envelope
attn = 70-intensity;
n = n.*10^(-attn/20); % attenuate
n = [n,zeros(1,(1/clickRate*fs)-length(n))];
n = repmat(n,1,nClicks);
stimArray(1,loc:loc+length(n)-1) = n;
stimArray(1,:) = conv(stimArray(1,:),FILT,'same');
events(1,loc:loc+length(n)-1) = ones(1,length(n))*5;
   
stimArray(2,:) = events;
stim = stimArray;

if any(stim>10)
    break
end
if any(stim<-10)
    break
end

stimInfo.seed = seed;
stimInfo.filename = filename;
stimInfo.fs = fs;
stimInfo.noiseDur = nDur; % noise duration in ms
stimInfo.nClicks = nClicks; % number of noise clicks
stimInfo.clickRate = clickRate; % in Hz
stimInfo.intensity = intensity; % in dB
stimInfo.ISI = ISI; % inter stim interval duration in ms
stimInfo.repeats = repeats; % number of repeats of each stim
stimInfo.filter = filterName;

chunk_size = []; nbits = 16;
fn = [filename '.wav'];
wavwrite_append((stim/10)', fn, chunk_size, fs, nbits)
save([filename '.mat'],'stimInfo')













