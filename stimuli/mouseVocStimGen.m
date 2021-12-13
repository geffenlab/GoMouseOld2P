% mouse vocalisation stim gen

clear all
seed = 42;
rng(seed);

filename = [datestr(now,'yyyymmdd') '_mouseVocs85dB_5reps' ];
fs = 400000;
vocID = [5, 6, 14, 18, 25, 34, 86, 101, 107, 71, 120]; %5=complex,6=chevron,8=complex,...
% 14=FMup,25=FMup,34=chevron, 86 = FMdown,101= FMdown,107=complex,
% 71&120=flat
ITI = 6; % inter tone interval duration in ms
attenuations = [-15];
vocrep = 10; % number of repeats of each voc per second
IVO = 1/vocrep; % inter-vocalisation onset in ms
repeats = 5; % number of reps


% Load the calibration filter
filtName = 'C:\calibration\Filters\20170420_2PspkrNidaqInvFilt_3k-80k_fs400k.mat';
load(filtName);

% Create stimuli
vocLoc = ['C:\Users\Chris\Dropbox\GeffenLab\Kath\mouseVocalisations\extracted\mats\'];
vocs = cell(length(vocID),1);
for ii=1:length(vocID)
    f = dir([vocLoc sprintf('mouseVoc%03d_*.mat',vocID(ii))]);
    load([vocLoc f.name]);
    vocs{ii} = (voc-mean(voc))/max(voc);
    vocs{ii} = vocs{ii}.*10^(-attenuations/20); % attenuate
    vocs{ii} = conv(vocs{ii},FILT,'same');
%     vocs{ii} = vocs{ii}(abs(vocs{ii})>2.4e-5);
end
    
stimArray = zeros(length(vocs),fs+(fs*6));
events = stimArray;
vHz = fs/vocrep;
for ii=1:length(vocs)
    for jj = 1:vocrep
        stimArray(ii,(jj-1)*vHz+1:(jj-1)*vHz+length(vocs{ii})) = vocs{ii};
        events(ii,(jj-1)*vHz+1:(jj-1)*vHz+length(vocs{ii})) = ones(1,length(vocs{ii}))*5;
    end
end

% for ii=1:length(vocs)
% %     spectrogram(stimArray(ii,:),200,20,1000,fs,'yaxis');
%       spectrogram(vocs{ii},200,20,1000,fs,'yaxis');
%     pause()
% end

% for ii=1:length(vocs)
%     plot(stimArray(ii,:))
%     hold on
%     plot(events(ii,:))
%     pause()
%     clf
% end
    

% Now create your vector
stim=zeros(repeats*size(stimArray,1)*size(stimArray,2),2);
ind = 1; % initiate index
order=[];
for ii = 1:repeats
    disp(ii)
    ro = randperm(size(stimArray,1),size(stimArray,1)); % select random order
    stimT = reshape(stimArray(ro,:)',1,length(stimArray(:)));
    evT = reshape(events(ro,:)',1,length(events(:)));
    stim((ii-1)*length(stimArray(:))+1:ii*length(stimArray(:)),1)=stimT;
    stim((ii-1)*length(events(:))+1:ii*length(events(:)),2)=evT;
    order = [order,ro];
end

% Make stim info
stimInfo.seed = seed;
stimInfo.filename = filename;
stimInfo.fs = fs;
stimInfo.vocIDs = vocID;
stimInfo.ISI = ITI; % inter stim interval duration in ms
stimInfo.vocrep = vocrep; % number of repeats of each voc per second
stimInfo.IVO = IVO; % inter-vocalisation onset in ms
stimInfo.attenuations = attenuations;
stimInfo.filterName = filtName;
stimInfo.index = vocID;
stimInfo.order = order;
stimInfo.stimGenFunc = 'mouseVocStimGen.m';
% order = toneOrder;

%%
chunk_size = []; nbits = 16;
fn = [filename '.wav'];
stim = (stim/10);
wavwrite_append(stim, fn, chunk_size, fs, nbits)
save([filename '_stimInfo.mat'],'stimInfo')
