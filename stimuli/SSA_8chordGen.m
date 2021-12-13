% Generate stimuli for FRA construction
% toneSequenceGen %(fs,minF,maxF,stepType,nSteps,octaveSteps,duration,ITI,attenuations,repeats)
clear all
clc

seed = 50;
rng(seed);

filename = ['E:\stimuli\Kath\SSA\' datestr(now,'yyyymmdd') 'harangTest'];
fs = 400000;
% fs = 40000;
f = [5000, 50000];
stepType = 'log'; % log or octaves, if length(f)>3 then those frequencies will be used
minF = f(1); % minimum frequency to test
maxF = f(2); % max frequency to test
octaveSteps = 0.5; % distance between the tones in octaves
nLogSteps = 20; % number of log steps
tDur = 150; % tone duration in ms
ICI = 650; % inter chord interval duration in ms
ISI = 4000; % inter stimulus interval
attenuation = 0; % 0 = 70 dB, -10 = 80dB, 10 = 60 dB
nTonesInChord = 10; % number of tones in the chords
totalChordsPerSeq = 8; % number of chords repeated



% Load the calibration filter
% filtName = 'C:\calibration\Filters\20170216_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
filtName = 'E:\calibration\Filters\20170711_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(filtName);

% Create stimuli
if length(f)==2
    if strcmp(stepType,'octaves')
        freqs = minF;
        while freqs(end)<maxF
            freqs(length(freqs)+1) = freqs(length(freqs))+(freqs(length(freqs))*octaveSteps);
        end
    elseif strcmp(stepType,'log')
        freqs = exp(linspace(log(minF),log(maxF),nLogSteps));
        freqs = round(freqs);
    end
else
    freqs = f;
end

% select two chords
r = freqs(randperm(length(freqs),nTonesInChord*2));
AB_f(1,:)= r(1:nTonesInChord);
AB_f(2,:) = r(nTonesInChord+1:end);

% Make the chords
ABtones = zeros(nTonesInChord,tDur/1000*fs,size(AB_f,1));
ABtones_sil = zeros(nTonesInChord,(tDur/1000*fs)+(ICI/1000*fs),size(AB_f,1));

silence = zeros(1,ICI/1000*fs); 
for jj=1:size(AB_f,1)
    for ii=1:nTonesInChord
            t = tone(AB_f(jj,ii),(3*pi)/2,tDur/1000,fs); % Make tone 
            t = t.*10^(-attenuation/20); % attenuate
            ABtones(ii,:,jj) = conv(t,FILT,'same');
            t = envelopeKCW(t,5,fs); % envelope envelope in ms
            ABtones_sil(ii,:,jj) = [ABtones(ii,:,jj),silence];
    end
end
clear ABtones
chords = squeeze(sum(ABtones_sil,1));
ev = ones(1,tDur/1000*fs)*5;
events = [ev,zeros(1,length(ABtones_sil(ii,:,jj))-length(ev))]';
% Make the sequences

AA = repmat(chords(:,1),8,1);
AB = [repmat(chords(:,1),7,1);chords(:,2)];
BB = repmat(chords(:,2),8,1);
BA = [repmat(chords(:,2),7,1);chords(:,1)];
events = repmat(events,8,1);

% preallocate variable
ISI_sil = zeros(ISI/1000*fs,1);
stim = [AA;ISI_sil;BB;ISI_sil;AB;ISI_sil;BA;ISI_sil];
stim(:,2) = repmat([events;ISI_sil],4,1);

% Make stim info
stimInfo.seed = seed;
stimInfo.filename = filename;
stimInfo.fs = fs;
stimInfo.frequencies = freqs;
stimInfo.stepType = stepType; % log or octaves, if length(f)>3 then those frequencies will be used
stimInfo.octaveSteps = octaveSteps; % distance between the tones in octaves
stimInfo.nLogSteps = nLogSteps; % number of log steps
stimInfo.tDur = tDur; % tone duration in ms
stimInfo.filterName = filtName;
stimInfo.nTonesInChord = nTonesInChord; % number of tones in the chords
stimInfo.ICI = ICI; % inter chord interval duration in ms
stimInfo.ISI = ISI; % inter stimulus interval
stimInfo.attenuation = attenuation; % 85dB
stimInfo.nTonesInChord = nTonesInChord; % number of tones in the chords
stimInfo.totalChordsPerSeq = totalChordsPerSeq; % number of chords repeated
stimInfo.index = {'AA','BB','AB','BA'};
stimInfo.chordFreqs = AB_f;


%%
stim = stim/10;
chunk_size = []; nbits = 16;
fn = [filename '.wav'];
wavwrite_append(stim, fn, chunk_size, fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')













