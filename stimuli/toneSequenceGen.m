% Generate stimuli for FRA construction
% toneSequenceGen %(fs,minF,maxF,stepType,nSteps,octaveSteps,duration,ITI,attenuations,repeats)
clear 
seed = round(rand(1)*1000);
rng(seed);
disp(seed)

filename = ['E:\stimuli\Kath\2P_FRA\' datestr(now,'yyyymmdd') '_2P_5k-32k_' sprintf('%03d',seed) ];
fs = 400000; % sample rate
% f = 2.^(linspace(log2(4000),log2(30000),3));
f = [5000,6145,7553,9283,11400,12492,13689,15000,17235,21184,26036,32000]; % Kath 2P stim
% CS = [11400 15000]; % if you do not want to add in conditioned tones, leave empty
stepType = 'log'; % log or octaves, if length(f)>3 then those frequencies will be used
minF = f(1); % minimum frequency to test
maxF = f(2); % max frequency to test
octaveSteps = 1/3; % distance between the tones in octaves
nLogSteps = 4; % number of log steps
tDur = 100; % tone duration in ms
ITI = 4000; % inter tone interval duration in ms
attenuations = [0];
repeats = 25; % number of repeats of EACH tone


% Load the calibration filter
filtName = 'E:\GitHub\filters\20191125_ABR_ephysBoothSpkr_300h-80k_fs400k.mat';
load(filtName);

% Create stimuli
totalDur = (ITI+tDur)/1000; % total duration in s
toneDur = tDur/1000;

if length(f)==2
    if strcmp(stepType,'octaves')
        freqs = minF;
        while freqs(end)<maxF
            freqs(length(freqs)+1) = freqs(length(freqs))+(freqs(length(freqs))*octaveSteps);
        end
        freqs = round(freqs(1:end-1));
    elseif strcmp(stepType,'log')
        freqs=exp(linspace(log(minF),log(maxF),nLogSteps));
        freqs = round(freqs);
    end
else
    freqs = f;
end

% preallocate variable
stimArray = zeros(length(freqs)*length(attenuations),round(totalDur*fs));
events = stimArray;
loc = 1; % placement of the tone in the zeros
ind = 1; % initiate index
for ii = 1:length(freqs)
    for jj = 1:length(attenuations)
        t = tone(freqs(ii),(3*pi)/2,toneDur,fs); % Make tone
        t = envelopeKCW(t,5,fs); % envelope % duration of envelope in ms
        t = t.*10^(-attenuations(jj)/20); % attenuate
        stimArray(ind,loc:loc+length(t)-1) = t;
%         stimArray(ind,:) = conv(stimArray(ind,:),FILT,'same');
        events(ind,loc:loc+length(t)-1) = ones(1,length(t))*5;
        index(ind,1) = freqs(ii);
        index(ind,2) = attenuations(jj);
        ind = ind+1;
    end
end

% Now create your vector
stim=zeros(repeats*size(stimArray,1)*size(stimArray,2),2);
ind = 1; % initiate index
toneOrder=[];
for ii = 1:repeats
    disp(ii)
    ro = randperm(size(stimArray,1),size(stimArray,1)); % select random order
    stimT = reshape(stimArray(ro,:)',1,length(stimArray(:)));
    evT = reshape(events(ro,:)',1,length(events(:)));
    stim((ii-1)*length(stimArray(:))+1:ii*length(stimArray(:)),1)=stimT;
    stim((ii-1)*length(events(:))+1:ii*length(events(:)),2)=evT;
    toneOrder = [toneOrder,ro];
end

stim(:,1) = conv(stim(:,1),FILT,'same'); % apply filter

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
stimInfo.attenuations = attenuations;
stimInfo.repeats = repeats; % number of repeats of each tone
stimInfo.filterName = filtName;
stimInfo.index = index;
stimInfo.order = toneOrder;
stimInfo.stimGenFunc = 'toneSequenceGen.m';
% order = toneOrder;
disp(round(length(stim)/fs))
disp([num2str(round(length(stim)/fs/60)) ' mins'])
%%
chunk_size = []; nbits = 16;
fn = [filename '.wav'];
stim = (stim/10);
wavwrite_append(stim, fn, chunk_size, fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')












