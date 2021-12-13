% Generate stimuli for FRA construction
% toneSequenceGen %(fs,minF,maxF,stepType,nSteps,octaveSteps,duration,ITI,attenuations,repeats)
clear all

filename = [datestr(now,'yyyymmdd') 'retigaCont_updownInterL_5-40k_90dB_4s' ];
fs = 400000;
tDur = 0.025; % tone duration in s
period = 4; % period duration
nTones = period/tDur;
maxF = 40000;
minF = 5000;
f = round(exp(linspace(log(minF),log(maxF),nTones)));
fup = round(exp(linspace(log(minF),log(maxF),nTones)));
fdown = round(exp(linspace(log(maxF),log(minF),nTones)));
ITI = 0; % inter tone interval duration in ms
ISI = 1; % interval between the up and down sweeps
attenuations = [-20];
repeats = 10; % number of repeats of each tone


% Load the calibration filter
filtName = 'C:\calibration\Filters\20170216_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(filtName);


% preallocate variable
totalDur = tDur;
stimArray = zeros(nTones*length(attenuations),round(totalDur*fs));
events = stimArray;
loc = 1; % placement of the tone in the zeros
ind = 1; % initiate index
eventDur = 0.01*fs;
for ii = 1:length(f)
    for jj = 1:length(attenuations)
        t = tone(f(ii),(3*pi)/2,tDur,fs); % Make tone
        t = t.*10^(-attenuations(jj)/20); % attenuate
        stimArray(ind,loc:loc+length(t)-1) = t;
        stimArray(ind,:) = envelopeKCW(stimArray(ind,:),5,fs); % envelope
%         stimArray(ind,:) = conv(stimArray(ind,:),FILT,'same');
        events(ind,loc:loc+eventDur-1) = ones(1,eventDur)*5;
        index(ind,1) = f(ii);
        index(ind,2) = attenuations(jj);
        ind = ind+1;
    end
end

% Now create your vector
% Do the up

stim=zeros(2,repeats*size(stimArray,1)*size(stimArray,2));
toneOrder=[];
for ii = 1:repeats
        disp(ii)
        ro=1:size(stimArray,1);
        stimT = reshape(stimArray',1,length(stimArray(:)));
        evT = reshape(events(ro,:)',1,length(events(:)));
        stim(1,(ii-1)*length(stimArray(:))+1:ii*length(stimArray(:)))=stimT;
        stim(2,(ii-1)*length(events(:))+1:ii*length(events(:)))=evT;
        toneOrder = [toneOrder,ro];
end
stim = [zeros(2,ISI*fs),stim];
% Do the down
stimArray = flipud(stimArray);
stim2=zeros(2,repeats*size(stimArray,1)*size(stimArray,2));
toneOrder2=[];
for ii = 1:repeats
        disp(ii)
        ro=1:size(stimArray,1);
        stimT = reshape(stimArray',1,length(stimArray(:)));
        evT = reshape(events(ro,:)',1,length(events(:)));
        stim2(1,(ii-1)*length(stimArray(:))+1:ii*length(stimArray(:)))=stimT;
        stim2(2,(ii-1)*length(events(:))+1:ii*length(events(:)))=evT;
        toneOrder2 = [toneOrder2,ro];
end
stim2 = [zeros(2,ISI*fs),stim2];
stim = [stim,stim2]'; clear stim2;
toneOrder = [toneOrder,toneOrder2];
stim(:,1) = conv(stim(:,1),FILT,'same');

% Make stim info
% stimInfo.seed = seed;
stimInfo.filename = filename;
stimInfo.fs = fs;
stimInfo.frequencies = f;
stimInfo.stepType = 'log'; % log or octaves, if length(f)>3 then those frequencies will be used
stimInfo.tDur = tDur; % tone duration in ms
stimInfo.ITI = ITI; % inter tone interval duration in ms
stimInfo.attenuations = attenuations;
stimInfo.sweepReps = repeats; % number of repeats of each tone
stimInfo.filterName = filtName;
stimInfo.index = index;
stimInfo.order = toneOrder;
stimInfo.stimGenFunc = 'contToneSequenceGen.m';
stimifno.period = period; % period duration
stimInfo.InterSweepInt = ISI; % interval between the up and down sweeps

%%
chunk_size = fs; nbits = 16;
fn = [filename '.wav'];
stim = (stim/10);
wavwrite_append(stim, fn, chunk_size, fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')













