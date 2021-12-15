%% set variables
clear all
delete(instrfindall)
global ps % this global is needed for the continuous nidaq presentation
repeats = 1;
bkg=5; % silence pre stim in s, make whole number

%% prepare the nidaq card
disp('Initialising Nidaq card')
fs=200000;
ps.s = daq.createSession('ni');
ps.s.Rate = fs;
ao1=addAnalogOutputChannel(ps.s,'dev1',0,'Voltage');
ao2=addAnalogOutputChannel(ps.s,'dev1',1,'Voltage');
ai1=addAnalogInputChannel(ps.s,'dev1',1,'Voltage');
ai2=addAnalogInputChannel(ps.s,'dev1',8,'Voltage');
ps.s.Channels(1).InputType = 'SingleEnded';
ps.s.Channels(2).InputType = 'SingleEnded';

ps.fs=ps.s.Rate; % check sample rate
if ps.fs~=fs
    disp('Sample rate not available')
    return
end

%% Set up thje sound card

sc.fs = 192e3;
InitializePsychSound(1);
pause(1); 
sc.h = PsychPortAudio('Open', [], 1, 3, sc.fs, 2);



%% prepare the stimulus
stim=audioread('quickTonotopy_5reps.wav');
stim=stim*10;
stim = repmat(stim,repeats,1); % number of repeats
wn = round((1-mod(length(stim)/fs,1))*fs);
if wn~=1
    stim = [stim;zeros(wn,2)];
end
triggerDuration = 500; % in samples
stim = [[zeros(bkg*fs,1),[ones(triggerDuration,1)*5;zeros((bkg*fs)-triggerDuration,1)]];stim];
bkgPadding = zeros(fs*3,2);
stim = [stim;bkgPadding];
ps.stim = zeros(fs,2,length(stim)/fs);
for ii=1:length(stim)/fs
    ps.stim(:,:,ii) = stim((ii-1)*fs+1:ii*fs,:);
end

% Read all sound files and create & fill one dynamic audiobuffer for
% each sound chunk:
buffer = [];
j = 0;
nfiles = size(ps.stim,3);
for i=1:nfiles
    audiodata = ps.stim(:,:,i);
    [samplecount, ninchannels] = size(audiodata);
    audiodata = repmat(transpose(audiodata), 2 / ninchannels, 1);
    buffer(end+1) = PsychPortAudio('CreateBuffer', [], audiodata); %#ok<AGROW>   
end

PsychPortAudio('UseSchedule', sc.h, 1);
for i=1:nfiles
    % Play buffer(i) from startSample 0.0 seconds to endSample 1.0 
    % seconds. Play one repetition of each soundbuffer...
    PsychPortAudio('AddToSchedule', sc.h, buffer(i), 1, 0.0, 1.0, 1);
end



ps.lh = addlistener(ps.s,'DataRequired',@presentStimCont);
ps.la = addlistener(ps.s,'DataAvailable',@acquireCont);
ps.data=zeros(length(stim),2);
ps.s.IsContinuous = true;
ps.counter=1;
ps.acqCounter=1;
disp(['stimulus duration: ' num2str((length(stim)-length(bkgPadding))/fs) ' s, ' num2str((length(stim)-length(bkgPadding))/fs/60) ' mins'])
disp(['29.6 fps : ' num2str(((length(stim)-length(bkgPadding))/fs)*29.6) ' frames']);
%% present the stimulus
ps.stim=zeros(size(ps.stim));
x = squeeze(ps.stim(:,:,1));
queueOutputData(ps.s,x);
PsychPortAudio('Start', sc.h, [], 0, 1);
ps.s.startBackground();  
% ps.data = [ps.data;ps.s.event.Data];