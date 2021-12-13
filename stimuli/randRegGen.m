% randRegTonePatternGenerator

clear

% for ff=1:2 % For reg-rand and rand-reg

% Set variables
seed = 5;
rng(seed);
fs = 400000; % number of samples per seconds
tonePipDur = 50; % ms
envDur = 5; % envelope duration (ms)
tpdfs = tonePipDur/1000*fs; % tone pip duration in samples
regCycle = 5; % number of tones that will repeat regularly
totalToneSeqs = 80;
totalTonePips = regCycle*totalToneSeqs; % total number of tone pips
nregseq = 3; % number of regular sequences
transtime = totalToneSeqs/2; % number of regular sequence repeats before transition
randtime = transtime-1; % number of random sequences repeats before transition
ISI = 0; % silence between stimuli
attenuation = -15;

% Load filter
filtName = 'C:\calibration\Filters\20170216_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(filtName);

% Random tone variables
nRand=15; % number of tones in the log space between the frequencies (could be equal or not to regCycle) - number of tones in pool
lowFreq=5000;%179; % low limit of frequencies to choose from
highFreq=50000;%7246; % high limit of frequencies to choose from
freqs=exp(linspace(log(lowFreq),log(highFreq),nRand));
freqs = round(freqs,-2); % tone pool

% Select regular sequences (A, B, C...)
seq=[];
for ii = 1:nregseq
    seq(ii,:) = (randperm(length(freqs),regCycle));
end

% Make regular sequence index for transitions
temp1 = 1:nregseq;
for ii = 1:nregseq
    temp2 = temp1(temp1~= ii);
    n = nregseq-1;
    seqIndex((ii-1)*n+1:ii*n,1) = ii;
    seqIndex((ii-1)*n+1:ii*n,2) = temp2;
end

% pre-make the tones for selection
tonematrix=[];
for ii = 1:length(freqs)
    t = tone(freqs(ii), 1, tonePipDur/1000, fs);
    t = envelopeKCW(t,envDur,fs);
    t = conv(t,FILT,'same');
    t = t.*10^(-attenuation/20); % attenuate
    tonematrix (ii,:) = t;
end


stim = zeros(size(seqIndex,1),totalTonePips*tpdfs);
stimP = [];
for ii = 1:size(seqIndex,1)
    
    for jj = 1:randtime
        r = randperm(length(freqs),regCycle);
        toneCycle=[];
        for kk = 1:regCycle
            toneCycle = [toneCycle,tonematrix(r(kk),:)];
        end
        stim(ii,((jj-1)*length(toneCycle)+1):jj*length(toneCycle)) = toneCycle;
        stimP(ii,(jj-1)*regCycle+1:jj*regCycle) = freqs(r);
    end
    
    toneCycle=[];
    for jj = 1:regCycle
        toneCycle = [toneCycle,tonematrix(seq(seqIndex(ii,1),jj),:)];
    end
    stim(ii,randtime*regCycle*tpdfs+1:(randtime+1)*regCycle*tpdfs) = toneCycle;
    stimP(ii,(randtime)*regCycle+1:(randtime+1)*regCycle) = freqs(seq(seqIndex(ii,2),:));
    
    % 10 consecutive regulars
    toneCycle=[];
    for jj = 1:regCycle
        toneCycle = [toneCycle,tonematrix(seq(seqIndex(ii,2),jj),:)];
    end
    stim(ii,(randtime+1)*regCycle*tpdfs+1:(totalTonePips*tpdfs)) = repmat(toneCycle,1,transtime);
    stimP(ii,(randtime+1)*regCycle+1:(transtime*2)*regCycle) = repmat(freqs(seq(seqIndex(ii,1),:)),1,transtime);
    
end

silence = zeros(1,ISI*fs);
x = [];  r = randperm(size(stimP,1),size(stimP,1));
% events = [];
events = repmat([ones(1,500),zeros(1,transtime*regCycle*tpdfs-500)],1,2*size(stim,1));
for ii = 1:size(stim,1)
    stim2 = [stim(r(ii),:), silence];
    x = [x,stim2];
%     events = [events,e];
end
stimP = stimP(r,:);
%
imagesc(stimP)
x = [x;events];
x = x/10; % /10 for saving as wav file

stimInfo.seed = seed;
stimInfo.fs = fs; % number of samples per seconds
stimInfo.tonePipDur = tonePipDur; % ms
stimInfo.envDur = envDur; % envelope duration (ms)
stimInfo.tpdfs = tpdfs; % tone pip duration in samples
stimInfo.regCycle = regCycle; % number of tones that will repeat regularly
stimInfo.totalToneSeqs = totalToneSeqs;
stimInfo.totalTonePips = totalTonePips; % total number of tone pips
stimInfo.nregseq = nregseq; % number of regular sequences
stimInfo.transtime = transtime; % number of regular sequence repeats before transition
stimInfo.randtime = randtime; % number of random sequences repeats before transition

% Random tone variables
stimInfo.nRand=nRand; % number of tones in the log space between the frequencies (could be equal or not to regCycle) - number of tones in pool
stimInfo.lowFreq=lowFreq;%179; % low limit of frequencies to choose from
stimInfo.highFreq=highFreq;%7246; % high limit of frequencies to choose from
stimInfo.freqs = freqs; % tone pool
stimInfo.ISI = ISI; % silence between stimuli in seconds
stimInfo.attenuation = attenuation;
stimInfo.filterName = filtName;
stimInfo.seqOrder = seqIndex(r,:);

% % sound(x/50,fs/4) % play the sound
% % spectrogram(stim(2,:), 1000, 0, 10000, fs,'yaxis');
fn = ['C:\Users\Chris\Documents\GitHub\Kath\stimuli\' datestr(now,'yyyymmdd') '_randReg_pool' num2str(nRand) '85dB_1rep.wav'];

%
chunk_size = []; nbits = 16;
wavwrite_append(x', fn, chunk_size, fs, nbits)
save(['C:\Users\Chris\Documents\GitHub\Kath\stimuli\' datestr(now,'yyyymmdd') '_randReg_pool' num2str(nRand) '85dB_1rep.mat'],'stimP','stimInfo')


% end



