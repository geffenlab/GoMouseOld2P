%% set variables
clear all
clear global
clc
delete(instrfindall)
global ps % this global is needed for the continuous nidaq presentation
repeats =5;
bkg=5; % silence pre stim in s, make whole number
stimFile = '20170308FRA_5-50kHz_4attns';
mouse = 'test';
ps.saveName = [datestr(now,'yymmdd') mouse '_' stimFile '_events.txt'];

cd('C:\Users\Chris\Documents\GitHub\Kath\2P_stim_PC\stimuli')

%% prepare the nidaq card
disp('Initialising Nidaq card')
fs=400000;
ps.s = daq.createSession('ni');
ps.s.Rate = fs;
ao1=addAnalogOutputChannel(ps.s,'dev1',0,'Voltage');
ao2=addAnalogOutputChannel(ps.s,'dev1',1,'Voltage');
ai1=addAnalogInputChannel(ps.s,'dev1',1,'Voltage');
ai2=addAnalogInputChannel(ps.s,'dev1',8,'Voltage');
ps.s.Channels(3).InputType = 'SingleEnded';
ps.s.Channels(4).InputType = 'SingleEnded';

ps.fs=ps.s.Rate; % check sample rate
if ps.fs~=fs
    disp('Sample rate not available')
    return
end
    
%% prepare the stimulus
stim=[];
stimInf = audioinfo([stimFile '.wav']);
if stimInf.SampleRate~=fs
    disp('STIM AT WRONG SAMPLE RATE!!')
    return
end
stimDur = stimInf.TotalSamples;
if stimDur*repeats/fs/60>30
    disp('STIM TOO LONG')
    keyboard
end

stim=audioread([stimFile '.wav']);

stim=stim*10;
% stim = repmat(stim,repeats,1); % number of repeats

% pad with zeros at the end if not a whole number of seconds
wn = fs-mod(length(stim)*repeats,fs); 
if wn~=fs
    zeroPad = zeros(wn,2);
else
    zeroPad=[];
end

% Add trigger to start recording
triggerDuration = 500; % in samples
triggerAcquisition = [[zeros(bkg*fs,1),[ones(triggerDuration,1)*5;zeros((bkg*fs)-triggerDuration,1)]]];

% Add some zeros at the end to make sure stimulus is not cut off
bkgPadding = zeros(fs*5,2);
% stim = [stim;bkgPadding];

% Divide the stim into 1 second chunks for continuous playback using the
% nidaq
ps.stim = zeros(fs,2,[length(stim)*repeats+length(zeroPad)+length(bkgPadding)+length(triggerAcquisition)]/fs);
% [triggerAcquisition; stim*repeats; bkgPadding]
stimIndex = 1;
for ii=1:length(triggerAcquisition)/fs
    ps.stim(:,:,stimIndex) = triggerAcquisition((ii-1)*fs+1:ii*fs,:);
    stimIndex = stimIndex+1;
end
clear triggerAcquisition
rm = [];
for jj=1:repeats
    if jj==repeats
        stim = [rm;stim;zeroPad];
    else
        stim = [rm;stim];
    end
    for ii=1:floor(length(stim)/fs)     
        ps.stim(:,:,stimIndex) = stim((ii-1)*fs+1:ii*fs,:);
        disp([num2str(ii) '/' num2str(length(stim)/fs)])
        stimIndex = stimIndex+1;
    end
    if mod(length(stim),fs)~=0
        rm2 = stim(end-mod(length(stim),fs)+1:end,:);
    else
        rm2=[];
    end
    stim = stim(length(rm)+1:end,:);
    rm=rm2;
end
clear stim

for ii=1:length(bkgPadding)/fs
    ps.stim(:,:,stimIndex) = bkgPadding((ii-1)*fs+1:ii*fs,:);
    stimIndex = stimIndex+1;
end
clear bkgPadding


% Add the continuous playback and acquisition listeners
sizeStim = size(ps.stim,1)*size(ps.stim,2)*size(ps.stim,3);
ps.lh = addlistener(ps.s,'DataRequired',@presentStimCont);
ps.la = addlistener(ps.s,'DataAvailable',@acquireCont_v2);
ps.fid = fopen(['C:\Experiments\newData\' ps.saveName],'W');
% ps.data=zeros(sizeStim/2,2); % preallocate recording data
ps.s.IsContinuous = true;
ps.counter=1; % start counters
ps.acqCounter=1; % start counters

% display stimulus info
disp(['stimulus duration: ' num2str(sizeStim/2/fs) ' s, ' num2str(sizeStim/2/fs/60) ' mins'])
fr = 29.534;
disp([num2str(fr) ' fps : ' num2str((sizeStim/2/fs)*fr) ' frames']); % estimate number of frames
%% present the stimulus
queueOutputData(ps.s,squeeze(ps.stim(:,:,1)));
ps.s.startBackground();  
% ps.data = [ps.data;ps.s.event.Data];
