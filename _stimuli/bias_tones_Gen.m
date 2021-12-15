%% Generate bias stimuli

clear 
s.seed = round(rand(1)*1000);
rng(s.seed);
disp(s.seed)

s.filename = ['E:\stimuli\Kath\bias\' datestr(now,'yyyymmdd') '_bias_' sprintf('%03d',s.seed) ];
s.fs = 400000; % sample rate
s.target_freq = 10000; % kHz
s.bias_octaves = [-1 -0.5 -0.25 -0.1 0.1 0.25 0.5 1]; % bias frequencies, octaves distance from s.target_freq
s.bias_freq = zeros(1,length(s.bias_octaves));
for ii = 1:length(s.bias_freq)
   s.bias_freq(ii) = 2.^(log2(s.target_freq) + s.bias_octaves(ii));
end
s.tDur = 100; % tone duration in ms
s.ITI = 300; % inter tone interval duration in ms
s.attenuations = [0];
s.bias_reps = 6;
s.repeats = 20; % number of s.repeats of EACH stim
s.ISI = 3000; % Duration between each stimulus


% Load the calibration filter
s.filtName = 'G:\GitHub\filters\20191108_2PspkrNidaqInvFilt_1-70k_fs400k.mat';
load(s.filtName);

% Create stimuli
totalDur = ((s.ITI+s.tDur)/1000)*(s.bias_reps+1); % total duration in s
toneDur = s.tDur/1000;

% preallocate variable
stimArray = zeros(length(s.bias_freq)*length(s.attenuations),round(totalDur*s.fs));
events = stimArray;

target = tone(s.target_freq,(3*pi)/2,toneDur,s.fs);
target = envelopeKCW(target,5,s.fs);
target = [target, zeros(1,(s.ITI/1000)*s.fs)];
ind = 1;
s.index = zeros(length(s.bias_freq)*length(s.attenuations),2);
for ii = 1:length(s.bias_freq)
    for jj = 1:length(s.attenuations)
        bias_tone = tone(s.bias_freq(ii),(3*pi)/2,toneDur,s.fs); % Make tone
        bias_tone = envelopeKCW(bias_tone,5,s.fs); % envelope % duration of envelope in ms
        bias_tone = repmat([bias_tone, zeros(1,(s.ITI/1000)*s.fs)],1, s.bias_reps);
        stim = [bias_tone, target];
        stim = stim.*10^(-s.attenuations(jj)/20); % attenuate
        stimArray(ind,:) = stim;
        events(ind,1:toneDur*s.fs) = ones(1,toneDur*s.fs)*5;
        s.index(ind,1) = s.bias_freq(ii);
        s.index(ind,2) = s.attenuations(jj);
        ind = ind+1;
    end
end

stimArray = [stimArray, zeros(size(stimArray,1),(s.ISI/1000)*s.fs)];
clear stim bias_tone target 

% Now create your vector
stim = zeros(s.repeats*size(stimArray,1)*size(stimArray,2),2);
ind = 1; % initiate s.index
s.stimIndexOrder = zeros(size(stimArray,1)*s.repeats,1);
for ii = 1:s.repeats
    disp(ii)
    ro = randperm(size(stimArray,1),size(stimArray,1)); % select random order
    stimT = reshape(stimArray(ro,:)',1,length(stimArray(:)));
    evT = reshape(events(ro,:)',1,length(events(:)));
    stim((ii-1)*length(stimArray(:))+1:ii*length(stimArray(:)),1) = stimT;
    stim((ii-1)*length(events(:))+1:ii*length(events(:)),2) = evT;
    s.stimIndexOrder((ii-1)*size(stimArray,1)+1:size(stimArray,1)*ii) = ro;
end

stim(:,1) = conv(stim(:,1),FILT,'same'); % apply filter

% Make stim info
s.stimGenFunc = 'bias_tones_Gen.m';
stimInfo = s;

disp([round(length(stim)/s.fs) ' seconds'])
disp([num2str((length(stim)/s.fs/60)) ' mins'])
%%
chunk_size = []; nbits = 16;
fn = [s.filename '.wav'];
stim = (stim/10);
wavwrite_append(stim, fn, chunk_size, s.fs, nbits)
save([s.filename '_stimInfo.mat'],'stimInfo')







