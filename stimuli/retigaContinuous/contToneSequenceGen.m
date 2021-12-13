% Generate stimuli for FRA construction
% toneSequenceGen %(fs,minF,maxF,stepType,nSteps,octaveSteps,duration,ITI,attenuations,repeats)
clear all


filename = [datestr(now,'yyyymmdd') 'retigaCont_40-5k_90dB_10s' ];
fs = 400000;
f = round(exp(linspace(log(40000),log(5000),400)));
stepType = 'log'; % log or octaves, if length(f)>3 then those frequencies will be used
minF = f(1); % minimum frequency to test
maxF = f(2); % max frequency to test
octaveSteps = 1/3; % distance between the tones in octaves
nLogSteps = 3; % number of log steps
tDur = 25; % tone duration in ms
ITI = 0; % inter tone interval duration in ms
attenuations = [-20];
repeats = 1; % number of repeats of each tone


% Load the calibration filter
filtName = 'C:\calibration\Filters\20170216_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
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
eventDur = 0.01*fs;
for ii = 1:length(freqs)
    for jj = 1:length(attenuations)
        t = tone(freqs(ii),(3*pi)/2,toneDur,fs); % Make tone
        t = t.*10^(-attenuations(jj)/20); % attenuate
        stimArray(ind,loc:loc+length(t)-1) = t;
        stimArray(ind,:) = envelopeKCW(stimArray(ind,:),5,fs); % envelope
%         stimArray(ind,:) = conv(stimArray(ind,:),FILT,'same');
        events(ind,loc:loc+eventDur-1) = ones(1,eventDur)*5;
        index(ind,1) = freqs(ii);
        index(ind,2) = attenuations(jj);
        ind = ind+1;
    end
end

% Now create your vector
% silenceVec = [zeros(1,preStimSil*fs);[ones(1,50),zeros(1,(preStimSil*fs)-50)]]*5;
stim=zeros(2,repeats*size(stimArray,1)*size(stimArray,2));
ind = 1; % initiate index
toneOrder=[];
for ii = 1:repeats
    disp(ii)
%     ro = randperm(size(stimArray,1),size(stimArray,1)); % select random order
     ro=1:size(stimArray,1);
    stimT = reshape(stimArray',1,length(stimArray(:)));
    evT = reshape(events(ro,:)',1,length(events(:)));
    stim(1,(ii-1)*length(stimArray(:))+1:ii*length(stimArray(:)))=stimT;
    stim(2,(ii-1)*length(events(:))+1:ii*length(events(:)))=evT;
    toneOrder = [toneOrder,ro];
end

stim(1,:) = conv(stim(1,:),FILT,'same');

% Make stim info
% stimInfo.seed = seed;
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
stimInfo.stimGenFunc = 'contToneSequenceGen.m';
% order = toneOrder;

%%
chunk_size = fs; nbits = 16;
fn = [filename '.wav'];
stim = (stim/10)';
wavwrite_append(stim, fn, chunk_size, fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')













