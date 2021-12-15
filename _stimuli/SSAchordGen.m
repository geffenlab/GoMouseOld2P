% Generate stimuli for FRA construction
% toneSequenceGen %(fs,minF,maxF,stepType,nSteps,octaveSteps,duration,ITI,attenuations,repeats)
clear all
clc

seed = 654;
rng(seed);

filename = ['E:\stimuli\Kath\SSA\' datestr(now,'yyyymmdd') 'ssaChords_A85B15_70dB_5-50k_' num2str(seed)];
fs = 400000;
f = [5000, 50000];
stepType = 'log'; % log or octaves, if length(f)>3 then those frequencies will be used
minF = f(1); % minimum frequency to test
maxF = f(2); % max frequency to test
octaveSteps = 0.5; % distance between the tones in octaves
nLogSteps = 20; % number of log steps
tDur = 150; % tone duration in ms
ISI = 850; % inter tone interval duration in ms
attenuations = [0];
repeats = NaN; % number of repeats of each tone
nTonesInChord = 7; % number of tones in the chords
percStandard = 85; % percentage of standards to oddballs
nChords = 100;


% Load the calibration filter
filtName = 'E:\calibration\Filters\20180426_2PspkrNidaqInvFilt_3-70k_fs400k.mat';
load(filtName);

% Create stimuli
if length(f)==2
    if strcmp(stepType,'octaves')
        freqs = minF;
        while freqs(end)<maxF
            freqs(length(freqs)+1) = freqs(length(freqs))+(freqs(length(freqs))*octaveSteps);
        end
    elseif strcmp(stepType,'log')
        freqs=exp(linspace(log(minF),log(maxF),nLogSteps));
        freqs = round(freqs);
    end
else
    freqs = f;
end

% select two chords
r = freqs(randperm(length(freqs),nTonesInChord*2));
AB(1,:)= r(1:nTonesInChord);
AB(2,:) = r(nTonesInChord+1:end);
oddballLocs = randperm(nChords,nChords*((100-percStandard)/100));
oddballLocs = sort(oddballLocs);

% Make the chords
ABtones = zeros(nTonesInChord,tDur/1000*fs,size(AB,1));
for jj=1:size(AB,1)
    for ii=1:nTonesInChord
        for kk=1:length(attenuations)
            t = tone(AB(jj,ii),(3*pi)/2,tDur/1000,fs); % Make tone
            t = envelopeKCW(t,5,fs); % envelope
            t = t.*10^(-attenuations(kk)/20); % attenuate
            ABtones(ii,:,jj) = conv(t,FILT,'same');
        end
    end
end

chords = squeeze(sum(ABtones,1));

% preallocate variable
stim = zeros(1,nChords*((tDur+ISI)/1000*fs)*length(attenuations));
events = stim;
loc = 1; % placement of the tone in the zeros
ind = 1; % initiate index
d = (ISI+tDur)/1000*fs;
for ii = 1:nChords
    if any(oddballLocs==ii)
        stim((ii-1)*d+1:(ii-1)*d+length(chords)) = chords(:,2);
    else
        stim((ii-1)*d+1:(ii-1)*d+length(chords)) = chords(:,1);
    end
    events((ii-1)*d+1:(ii-1)*d+length(chords)) = ones(1,length(chords))*5;
end

% Make stim info
stimInfo.seed = seed;
stimInfo.filename = filename;
stimInfo.fs = fs;
stimInfo.frequencies = freqs;
stimInfo.stepType = stepType; % log or octaves, if length(f)>3 then those frequencies will be used
stimInfo.octaveSteps = octaveSteps; % distance between the tones in octaves
stimInfo.nLogSteps = nLogSteps; % number of log steps
stimInfo.tDur = tDur; % tone duration in ms
stimInfo.ISI = ISI; % inter tone interval duration in ms
stimInfo.attenuations = attenuations;
stimInfo.repeats = repeats; % number of repeats of each tone
stimInfo.filterName = filtName;
stimInfo.nTonesInChord = 7; % number of tones in the chords
stimInfo.percStandard = 90; % percentage of standards to oddballs
stimInfo.nChords = 200;
stimInfo.index = AB;
order = ones(nChords,1);
order(oddballLocs) = 2;
stimInfo.order = order;

%%
stim = [stim;events];
stim = stim/10;
chunk_size = []; nbits = 16;
fn = [filename '.wav'];
wavwrite_append(stim', fn, chunk_size, fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')













