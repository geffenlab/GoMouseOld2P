% Generate stimuli for rate level function

clear
seed = 21;
rng(seed);

filename = 'E:\stimuli\Kath\widefield_noise\20171207_noiseBursts_50-90dBx20_400k';
fs = 400000;
nDur = 100; % noise duration in ms
ISI = 3000; % inter stim interval duration in ms
attenuations = [-20:20:20];
repeats = 20; % number of repeats of each tone

% Load the calibration filter
filterName = 'E:\calibration\Filters\20170711_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(filterName);

% Create stimuli
totalDur = (ISI+nDur)/1000; % total duration in s

% preallocate variable
stimArray = zeros(length(attenuations),round(totalDur*fs));
events = stimArray;
loc = 1; % placement of the tone in the zeros
ind = 1; % initiate index
for jj = 1:length(attenuations)
    n = rand(1,(nDur/1000)*fs);
    n = envelopeKCW(n,5,fs); % envelope
    n = n.*10^(-attenuations(jj)/20); % attenuate
    stimArray(ind,loc:loc+length(n)-1) = n;
    stimArray(ind,:) = conv(stimArray(ind,:),FILT,'same');
    events(ind,loc:loc+length(n)-1) = ones(1,length(n))*5;
    index(ind,1) = attenuations(jj);
    ind = ind+1;
end


% Now create your vector
% silenceVec = [zeros(1,preStimSil*fs);[ones(1,50),zeros(1,(preStimSil*fs)-50)]]*5;
stim=zeros(2,repeats*length(stimArray(:)));
ind = 1; % initiate index
attnOrder=[];
for ii = 1:repeats
    disp(ii)
    ro = randperm(size(stimArray,1),size(stimArray,1)); % select random order
    stimT = reshape(stimArray(ro,:)',1,length(stimArray(:)));
    evT = reshape(events(ro,:)',1,length(events(:)));
    stim(1,(ii-1)*length(stimArray(:))+1:ii*length(stimArray(:)))=stimT;
    stim(2,(ii-1)*length(events(:))+1:ii*length(events(:)))=evT;
    attnOrder = [attnOrder,ro];
end

stimInfo.seed = seed;
stimInfo.fs = fs;
stimInfo.noiseDur = nDur; % noise duration in ms
stimInfo.ISI = ISI; % inter stim interval duration in ms
stimInfo.attenuations = attenuations;
stimInfo.repeats = repeats; % number of repeats of each tone
stimInfo.filterName = filterName;
stimInfo.order = attnOrder;
stimInfo.index = index;
stimInfo.stimGenFunc = 'rateLevelNoiseGen.m';

disp(length(stim)/fs/60);
%%
chunk_size = []; nbits = 16;
fn = [filename '.wav'];
wavwrite_append((stim/10)', fn, chunk_size, fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')













