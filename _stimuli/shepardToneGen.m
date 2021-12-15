% Generate stimuli for FRA construction
% toneSequenceGen %(fs,minF,maxF,stepType,nSteps,octaveSteps,duration,ITI,attenuations,repeats)
clear all
seed = round(rand(1)*1000);
rng(seed);
disp(seed)

filename = ['E:\stimuli\Kath\shepard\' datestr(now,'yyyymmdd') '_shepard_70dB_' sprintf('%03d',seed) ];
fs = 400000; % sample rate

chordFreqs = [10000 22500];
% CS = [11400 15000]; % if you do not want to add in conditioned tones, leave empty
stepType = 'log'; % log or octaves, if length(f)>3 then those frequencies will be used
% minF = f(1); % minimum frequency to test
% maxF = f(2); % max frequency to test
octaveSteps = 1/2; % distance between the tones in octaves
nLogSteps = 4; % number of log steps
tDur = 0.3; % tone duration in ms
ITI = 1;
ISI = 5; % inter tone interval duration in ms
attenuations = [0];
repeats = 20; % number of repeats of each tone
nContextTones = 10;


% Load the calibration filter
filtName = 'E:\calibration\Filters\20170711_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
% filtName = 'E:\calibration\Filters\20160825_2PspkrNidaqInvFilt_3k-70k_fs200k.mat';
load(filtName);
ci = [1 2; 2 3; 1 3]; % low, high, equal context
ind = 1;
stim = zeros((length(chordFreqs)-1)*3,(ITI+tDur)*fs*nContextTones + (ITI+tDur)*2*fs + ISI*fs);
events = stim; ev = ones(1,0.05*fs)*5;
for ii = 1:length(chordFreqs)-1
    f = log2space(chordFreqs(1),chordFreqs(2),3);
    for jj = 1:3
        events(ind,1:length(ev)) = ev;
        cf = log2space(f(ci(jj,1)),f(ci(jj,2)),10);
        rt = cf(randperm(length(cf)));
        for tt = 1:length(rt)
            t = tone(rt(tt),1,tDur,fs);
            t = conv(t,FILT,'same');
            stim(ind,(tt-1)*(ITI+tDur)*fs+1:((tt-1)*(ITI+tDur)*fs)+fs*tDur) = envelopeKCW(t,5,fs);
        end
        tt = tt+1;
        t = tone(f(2),1,tDur,fs);
        t = conv(t,FILT,'same');
        stim(ind,(tt-1)*(ITI+tDur)*fs+1:((tt-1)*(ITI+tDur)*fs)+fs*tDur) = envelopeKCW(t,5,fs);
        tt = tt+1;
        t = (tone(f(1),1,tDur,fs) + tone(f(3),1,tDur,fs));
        t = t/max(abs(t));
        t = conv(t,FILT,'same');
        stim(ind,(tt-1)*(ITI+tDur)*fs+1:((tt-1)*(ITI+tDur)*fs)+fs*tDur) = envelopeKCW(t,5,fs);
        ind = ind+1;
    end
end
   
stimRep = zeros(length(stim(:))*repeats,2);
order = [];
for ii = 1:repeats
    r = randperm(size(stim,1));
    order = [order,r];
    stimRep((ii-1)*length(stim(:))+1:ii*length(stim(:)),1) = reshape(stim(r,:)',1,[]);
    stimRep((ii-1)*length(stim(:))+1:ii*length(stim(:)),2) = reshape(events(r,:)',1,[]);
end
    

% Make stim info
stimInfo.seed = seed;
stimInfo.filename = filename;
stimInfo.fs = fs;
stimInfo.frequencies = f;
stimInfo.stepType = stepType; % log or octaves, if length(f)>3 then those frequencies will be used
stimInfo.octaveSteps = octaveSteps; % distance between the tones in octaves
stimInfo.nLogSteps = nLogSteps; % number of log steps
stimInfo.tDur = tDur; % tone duration in ms
stimInfo.ITI = ITI; % inter tone interval duration in ms
stimInfo.ISI = ISI; % inter stimulus interval duration in ms
stimInfo.attenuations = attenuations;
stimInfo.repeats = repeats; % number of repeats of each tone
stimInfo.filterName = filtName;
% stimInfo.index = index;
stimInfo.order = order;
stimInfo.stimGenFunc = 'shepardToneGen.m';
stimInfo.chordFreqs = chordFreqs;
stimInfo.nContextTones = nContextTones;
stimInfo.contexts = {'low','high','equal'};

disp(round(length(stimRep)/fs))
disp([num2str(round(length(stimRep)/fs/60)) ' mins'])
%%
chunk_size = []; nbits = 16;
fn = [filename '.wav'];
stim = (stimRep/10);
wavwrite_append(stim, fn, chunk_size, fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')













