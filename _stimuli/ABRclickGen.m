clear

fs = 400000;
% click = -1:1/(fs/1000):1;
clickDur = fs/1000; % 1 ms
click = rand(1,clickDur);
ITI = (0.05*fs)-length(click);

stim = [click zeros(1,ITI)]';
stim(:,2) = [ones(1,0.01*fs)*5 zeros(1,length(stim)-(0.01*fs))]';
% stim(:,1) = stim(:,1).*10^(10/20); % 80 dB
 
% Load the calibration filter
filtName = 'E:\GitHub\filters\20191125_ABR_ephysBoothSpkr_300h-80k_fs400k.mat';
% filtName = 'E:\calibration\Filters\20160825_2PspkrNidaqInvFilt_3k-70k_fs200k.mat';
load(filtName);
stim(:,1) = conv(stim(:,1),FILT,'same'); % apply filter
stim = repmat(stim,40,1);
fn = ['H:\ABR\20201008_ABRclick_40reps_1ms'];
stim = (stim/10);
audiowrite([fn '.wav'],stim,fs,'BitsPerSample',64);
stimInfo.fs = fs;
stimInfo.filterName = filtName;
save([fn '_stimInfo.mat'],'stimInfo')

%% new code:
clear
rng(9)
fs = 400000;
clickDur = fs/1000; % 1 ms
clickRate = 20; % Hz
reps = 40;
totalDur = (fs/clickRate)*reps;
ITI = (0.05*fs)-clickDur;
singleClickDur = clickDur+ITI;
event = ones(0.01*fs,1)*5;
stim = zeros(totalDur,2);

for ii = 1:reps
    click = rand(clickDur,1);
    ind(1) = (ii-1)*singleClickDur+1;
    ind(2) = (ii-1)*singleClickDur+clickDur;
    ind(3) = (ii-1)*singleClickDur+length(event);
    stim(ind(1):ind(2),1) = click;
    stim(ind(1):ind(3),2) = event;
end

% Load the calibration filter
filtName = 'E:\GitHub\filters\20191125_ABR_ephysBoothSpkr_300h-80k_fs400k.mat';
load(filtName);
stim(:,1) = conv(stim(:,1),FILT,'same'); % apply filter

%% SAVE
fn = ['H:\ABR\20201008_ABRclick_40reps_1ms'];
stim = (stim/10);
audiowrite([fn '.wav'],stim,fs,'BitsPerSample',64);
stimInfo.fs = fs;
stimInfo.filterName = filtName;
save([fn '_stimInfo.mat'],'stimInfo')


