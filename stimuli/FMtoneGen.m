% Generate stimuli for FM tones

clear
seed = 183;
rng(seed);
addpath('C:\Users\Chris\Documents\GitHub\Kath\FM_stimuli\');
filename = [datestr(now,'yyyymmdd') 'FMtones_90dB_15kcarrier5kmod10-80modRate_400k'];

% Variables
fs = 400000;
carrierFreq = 15000; % Carrier Frequency
phi = 3*pi/2; % phase of the carrier frequency in radians, can change this to vary phase that the stimulus starts from
fModRate = [0:10:100]; % modulation rate
fModDepth = 5000; % depth of the modulation in Hz
AM = 0; % change to 1 if you also want AM moduation
tDur = 300; % tone duration in ms
ITI = 5000; % inter tone interval duration in ms
attenuations = [-20];
repeats = 1; % number of repeats of each tone

% Load the calibration filter
filtName = 'C:\calibration\Filters\20170216_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(filtName);

% Create stimuli
totalDur = (ITI+tDur)/1000; % total duration in s
toneDur = tDur/1000;

% preallocate variable
stimArray = zeros(length(fModRate)*length(attenuations),round(totalDur*fs));
events = stimArray;
loc = 1; % placement of the tone in the zeros
ind = 1; % initiate index
for ii = 1:length(fModRate)
    for jj = 1:length(attenuations)
        t = FM_stimGen(fs,carrierFreq,fModRate(ii),toneDur,fModDepth,phi); % Make tone
        t = envelopeKCW(t,10,fs); % envelope
        t = t.*10^(-attenuations(jj)/20); % attenuate
        stimArray(ind,loc:loc+length(t)-1) = t;
        stimArray(ind,:) = conv(stimArray(ind,:),FILT,'same');
        events(ind,loc:loc+length(t)-1) = ones(1,length(t))*5;
        index(ind,1) = fModRate(ii);
        index(ind,2) = attenuations(jj);
        ind = ind+1;
    end
end

% Now create your vector
% silenceVec = [zeros(1,preStimSil*fs);[ones(1,50),zeros(1,(preStimSil*fs)-50)]]*5;
stim=zeros(2,repeats*length(stimArray(:)));
ind = 1; % initiate index
order=[];
for ii = 1:repeats
    disp(ii)
    ro = randperm(size(stimArray,1),size(stimArray,1)); % select random order
    stimT = reshape(stimArray(ro,:)',1,length(stimArray(:)));
    evT = reshape(events(ro,:)',1,length(events(:)));
    stim(1,(ii-1)*length(stimArray(:))+1:ii*length(stimArray(:)))=stimT;
    stim(2,(ii-1)*length(events(:))+1:ii*length(events(:)))=evT;
    order = [order,ro];
end

% Make stim info
stimInfo.seed = seed;
stimInfo.filename = filename;
stimInfo.fs = fs;
stimInfo.carrierFreq = carrierFreq; % Carrier Frequency
stimInfo.phi = phi; % phase of the carrier frequency in radians, can change this to vary phase that the stimulus starts from
stimInfo.fModRate = fModRate; % modulation rate
stimInfo.fModDepth = fModDepth; % depth of the modulation in Hz
stimInfo.AM = AM; % change to 1 if you also want AM moduation
stimInfo.stimDur = tDur; % tone duration in ms
stimInfo.ISI = ITI; % inter tone interval duration in ms
stimInfo.attenuations = attenuations;
stimInfo.repeats = repeats; % number of repeats of each tone
stimInfo.filterName = filtName;
stimInfo.index = index;
stimInfo.order = order;

chunk_size = []; nbits = 16;
fn = [filename '.wav'];
wavwrite_append((stim/10)', fn, chunk_size, fs, nbits)
save([filename '.mat'],'stimInfo')













