% Generate stimuli for FRA construction
% toneSequenceGen %(fs,minF,maxF,stepType,nSteps,octaveSteps,duration,ITI,attenuations,repeats)
clear all
seed = round(rand(1)*1000);
rng(seed);
disp(seed)
stimInfo.seed = seed;

filename = ['E:\stimuli\Kath\widefield_tones\' datestr(now,'yyyymmdd') '_SAMTonesX10_300ms_5-30k_70-50dB_' sprintf('%03d',seed) ];
stimInfo.fs = 400000;
stimInfo.f = [5000, 15000, 30000];
stimInfo.stepType = 'log'; % log or octaves, if length(f)>3 then those frequencies will be used
stimInfo.octaveSteps = 0.5; % distance between the tones in octaves
stimInfo.nLogSteps = 5; % number of log steps
stimInfo.tDur = 300; % tone duration in ms
stimInfo.ITI = 3000; % inter tone interval duration in ms
stimInfo.attenuations = [0:10:20]; % dB
stimInfo.repeats = 10; % number of repeats of each tone
stimInfo.SAMrate = [10]; % Hz


% Load the calibration filter
stimInfo.filtName = 'E:\calibration\Filters\20170711_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(stimInfo.filtName);

% Create stimuli
totalDur = (stimInfo.ITI+stimInfo.tDur)/1000; % total duration in s
toneDur = stimInfo.tDur/1000;
minF = stimInfo.f(1); % minimum frequency to test
maxF = stimInfo.f(2); % max frequency to test
if length(stimInfo.f)==2
    if strcmp(stimInfo.stepType,'octaves')
        freqs = minF;
        while freqs(end)<maxF
            freqs(length(freqs)+1) = freqs(length(freqs))+(freqs(length(freqs))*stimInfo.octaveSteps);
        end
    elseif strcmp(stimInfo.stepType,'log')
        freqs=exp(linspace(log(minF),log(maxF),stimInfo.nLogSteps));
        freqs = round(freqs);
    end
else
    freqs = stimInfo.f;
end

% Make sam 
for ii = 1:length(stimInfo.SAMrate)
    sam(ii,:) = tone(stimInfo.SAMrate(ii)/2,-pi/2,toneDur,stimInfo.fs); % Make tone
end
 
% preallocate variable
stimArray = zeros(length(freqs)*length(stimInfo.attenuations)*length(stimInfo.SAMrate),round(totalDur*stimInfo.fs));
events = stimArray;
loc = stimInfo.ITI/2000*stimInfo.fs; % placement of the tone in the zeros
ind = 1; % initiate index
for ii = 1:length(freqs)
    for jj = 1:length(stimInfo.attenuations)
        for kk=1:length(stimInfo.SAMrate)
            t = tone(freqs(ii),(3*pi)/2,toneDur,stimInfo.fs); % Make tone
            t = t.*sam(kk,:);
            t = envelopeKCW(t,10,stimInfo.fs); % envelope
            t = t.*10^(-stimInfo.attenuations(jj)/20); % attenuate
            stimArray(ind,loc:loc+length(t)-1) = t;
%             stimArray(ind,:) = conv(stimArray(ind,:),FILT,'same');
            events(ind,loc:loc+length(t)-1) = ones(1,length(t))*5;
            index(ind,1) = freqs(ii);
            index(ind,2) = stimInfo.attenuations(jj);
            index(ind,3) = stimInfo.SAMrate(kk);
            ind = ind+1;
        end
    end
end

% Now create your vector
% silenceVec = [zeros(1,preStimSil*fs);[ones(1,50),zeros(1,(preStimSil*fs)-50)]]*5;
stim=zeros(stimInfo.repeats*length(stimArray(:)),2);
ind = 1; % initiate index
toneOrder=[];
for ii = 1:stimInfo.repeats
    disp(ii)
    ro = randperm(size(stimArray,1),size(stimArray,1)); % select random order
    stimT = reshape(stimArray(ro,:)',1,length(stimArray(:)));
    evT = reshape(events(ro,:)',1,length(events(:)));
    stim((ii-1)*length(stimArray(:))+1:ii*length(stimArray(:)),1)=stimT;
    stim((ii-1)*length(events(:))+1:ii*length(events(:)),2)=evT;
    toneOrder = [toneOrder,ro];
end

stim(:,1) = conv(stim(:,1),FILT,'same');

stimInfo.index = index;
stimInfo.order = toneOrder;
stimInfo.stimGenFunc = 'SAMtoneSequenceGen.m';

disp(['Stim duration: ' num2str(length(stim)/stimInfo.fs) ' s'])
disp(['Stim duration: ' num2str(length(stim)/stimInfo.fs/60) ' mins'])
%%
chunk_size = []; nbits = 16;
fn = [filename '.wav'];
stim = (stim/10);
wavwrite_append(stim, fn, chunk_size, stimInfo.fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')













