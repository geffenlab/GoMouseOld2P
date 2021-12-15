% Generate stimuli for rate level function

clear
seed = 352;
rng(seed);

filename = 'E:\stimuli\Kath\narrowbandNoise\20171205_widefield';
fs = 400000;
nDur = 0.25; % noise duration in ms
ISI = 2; % inter stim interval duration in ms
bandwidths = [3.75 7.5; 7.5 15; 15 30; 30 60]; % noise bandwidths
attenuations = [-10 5 20];
repeats = 10; % number of repeats of each noise busrt

% Load the calibration filter
filterName = 'E:\calibration\Filters\20170711_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(filterName);

% Create stimuli
totalDur = (ISI+nDur)*repeats*length(attenuations)*size(bandwidths,1); % total duration in s

% Make the filters
for ii = 1:size(bandwidths,1)
    L = bandwidths(ii,1)*1000/(fs/2);
    H = bandwidths(ii,2)*1000/(fs/2);
    [b(ii,:),a(ii,:)] = butter(2,[L H]);
end

% preallocate variables
stim = zeros((nDur+ISI)*fs*repeats*size(bandwidths,1)*length(attenuations),2);
ev = [ones(nDur*fs,1)*5;zeros(ISI*fs,1)];
locs  = [(((1:repeats*size(bandwidths,1)*length(attenuations))*((ISI+nDur)*fs))-(((ISI+nDur)*fs)-1))',...
    ((1:repeats*size(bandwidths,1)*length(attenuations))*((ISI+nDur)*fs))'];
ind = 1; index = [];
for kk = 1:repeats
    disp(kk);
    stimArray = zeros(size(bandwidths,1),(nDur+ISI)*fs,length(attenuations));
    aord = []; bord=[];
    for jj = 1:length(attenuations)
        aord = [aord;ones(size(bandwidths,1),1)*jj];
        bord = [bord;(1:size(bandwidths,1))'];
        for ii = 1:size(bandwidths,1)
            n = rand(1,nDur*fs);
            n = n+1; % normalise mean to 1
            n = n.*10^(-attenuations(jj)/20); % attenuate
            n = filter(b(ii,:),a(ii,:),n); % band pass filter
            n = conv(n,FILT,'same'); % speaker filter
            n = envelopeKCW(n,5,fs); % envelope
            stimArray(ii,:,jj) = [n,zeros(1,ISI*fs)];
        end
    end
    ord = [bord,aord];
    ro = randperm(size(ord,1))';
    ord = ord(ro,:);
    for ii = 1:size(ord,1)
        stim(locs(ind,1):locs(ind,2),:) = [(stimArray(ord(ii,1),:,ord(ii,2)))',ev];
        ind = ind+1;
    end
    index = [index; ro];
end

order = [bord,aord];



stimInfo.seed = seed;
stimInfo.fs = fs;
stimInfo.noiseDur = nDur; % noise duration in s
stimInfo.ISI = ISI; % inter stim interval duration in s
stimInfo.bandwidths = bandwidths;
stimInfo.attenuations = attenuations;
stimInfo.repeats = repeats; % number of repeats of each noise
stimInfo.filterName = filterName;
stimInfo.order = order;
stimInfo.orderInfo = {'bandwidths','attenuations'};
stimInfo.index = index;
stimInfo.stimGenFunc = 'bandPassNoiseStimGen.m';

disp(length(stim)/fs/60);
%%
chunk_size = []; nbits = 16;
fn = [filename '.wav'];
wavwrite_append((stim/10), fn, chunk_size, fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')













