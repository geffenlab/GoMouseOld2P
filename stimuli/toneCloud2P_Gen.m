% Make tone clouds
clear all

seed = round(rand(1)*1000);
rng(seed);
disp(seed)

filename = ['E:\stimuli\Kath\toneClouds\' datestr(now,'yyyymmdd') '_toneClouds_1s_5k-32k_65dB_' sprintf('%03d',seed) ];
% set variables
stimInfo.fs              = 4e5;              % sample rate
stimInfo.tonePipDur      = 0.030;            % duration of each tone pip in cloud
stimInfo.totalDur        = 1;                % total duration of each tone cloud
L2S = log2space(5000,32000,8); x = zeros(length(L2S)-1,2);
for ii = 1:length(L2S)-1
    x(ii,:) = [L2S(ii) L2S(ii+1)];
end
stimInfo.cloudRange = [x; 5000 32000];       % tone cloud bandwidths
stimInfo.nLogSteps       = 10;               % number of tones in the range
stimInfo.envDur          = 0.005;            % duration of tone pip envelope
stimInfo.tonePipRate     = 100;              % presentation rate in Hz (determines tone overlap)
stimInfo.toneLevel       = 65;               % levels of tones in dB
stimInfo.repeats         = 20;
stimInfo.ISI             = 4;
stimInfo.nTonesPerCloud = stimInfo.totalDur*stimInfo.tonePipRate; % number of tones in stimulus
atten = 70-stimInfo.toneLevel; % convert tone level to attenuation from 70 dB (filter is set to make sounds at 70 dB)

cloudTones = zeros(size(stimInfo.cloudRange,1),stimInfo.nLogSteps);
for cc = 1:size(stimInfo.cloudRange,1)
    cloudTones(cc,:) = log2space(stimInfo.cloudRange(cc,1),stimInfo.cloudRange(cc,2),stimInfo.nLogSteps); % tones frequencies
end
  
% make matrix of tones
toneMat = zeros(size(cloudTones,1),size(cloudTones,2),stimInfo.tonePipDur*stimInfo.fs); % preallocate tone matrix
for ii = 1:size(cloudTones,1)
    for jj = 1:size(cloudTones,2)
        t = tone(cloudTones(ii,jj),3/2*pi,stimInfo.tonePipDur,stimInfo.fs); % make tone
        toneMat(ii,jj,:)  = envelopeKCW(t,stimInfo.envDur*1000,stimInfo.fs); % envelope it
    end
end


% make tone cloud
stimInfo.toneOrder = zeros(stimInfo.nTonesPerCloud,size(stimInfo.cloudRange,1),stimInfo.repeats);
stimInfo.rangeOrder = zeros(size(stimInfo.cloudRange,1),stimInfo.repeats);
oneRepDur = size(stimInfo.cloudRange,1)*(stimInfo.totalDur+stimInfo.ISI)*stimInfo.fs;
stim = zeros(oneRepDur*stimInfo.repeats,1);
events = stim; TTL = ones(round(0.025*stimInfo.fs),1)*5;
for ii = 1:stimInfo.repeats
    r2 = randperm(size(stimInfo.cloudRange,1));
    stimInfo.rangeOrder(:,ii) = r2;
    for jj = 1:size(stimInfo.cloudRange,1)
        r = randi(size(cloudTones,2),stimInfo.nTonesPerCloud,1); % randomly select tones from the range
        stimInfo.toneOrder(:,jj,ii) = r; % save order of tones ?
        tonePipOnsets = ((jj-1)*(stimInfo.totalDur+stimInfo.ISI)*stimInfo.fs) + (1:round(stimInfo.totalDur*stimInfo.fs)/stimInfo.tonePipRate:round(stimInfo.totalDur*stimInfo.fs)); % tone pip onsets
        tonePipOnsets = tonePipOnsets + ((ii-1)*oneRepDur);
        for kk = 1:stimInfo.nTonesPerCloud
            stim(tonePipOnsets(kk):tonePipOnsets(kk)+(stimInfo.tonePipDur*stimInfo.fs)-1,1)...
                = stim(tonePipOnsets(kk):tonePipOnsets(kk)+(stimInfo.tonePipDur*stimInfo.fs)-1,1) + squeeze(toneMat(r2(jj),r(kk),:));
        end
        events(tonePipOnsets(1):tonePipOnsets(1)+length(TTL)-1) = TTL;
    end
end

% stim(1:round(0.1*stimInfo.fs),2) = ones(round(0.1*stimInfo.fs),1)*5; % event time
stim  = stim.*10^(-atten/20); % attenuate   
stimInfo.filtName = 'E:\calibration\Filters\20170711_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(stimInfo.filtName);
stim = conv(stim,FILT,'same');

stim = [stim,events];

% check stim
% spectrogram(stim(:,1), 1000, 0, 10000, stimInfo.fs,'yaxis');

stimInfo.stimGenFunc = 'toneCloud2P_Gen.m';

disp(length(stim)/stimInfo.fs/60)

%% save stim
chunk_size = []; nbits = 16;
fn = [filename '.wav'];
stim = (stim/10);
wavwrite_append(stim, fn, chunk_size, stimInfo.fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')



    
    