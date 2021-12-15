%% present stimuli
cd('E:\GitHub\Kath\2P_stim_PC\')


delete(instrfindall)
daqreset
fclose('all');
clear all
clear global
clear
clc
global ps % this global is needed for the continuous nidaq presentation

% SET VARIABLES
mouse = 'K057';
ps.stimFiles(1).name = ['20170414_silence_60s' '.wav']; % files that you want to present (or comment if want to present many files contained in a folder (see below)
repeats = 20; 
ps.yaw = 350; % 2P microscope values
ps.pitch = 45; % 2P microscope values

% any notes
% exptInfo.notes = 'This objective angle is still not right and the mouse is moving/struggling a lot';


bkg=10; % silence pre stim in s, make whole number
% ps.stimFolder = 'E:\stimuli\Kath\current\fear\';
ps.stimFolder = 'E:\stimuli\Kath\current\fear\';
% ps.stimFiles = dir([ps.stimFolder '*.wav']); % uncomment if you want to
% run multiple stimuli contained in a folder
for ii=2:repeats
   ps.stimFiles(ii).name =  ps.stimFiles(1).name;
end

ps.saveName = [datestr(now,'yymmdd_HHMM') mouse '_' ps.stimFiles(1).name(1:end-4) '.txt']; % name for saving the acquired data (frame and event times)
ps.nFiles  = length(ps.stimFiles);
ps.mouse = mouse;
ps.reps = repeats;
ps.preStimSil = bkg;

%% prepare the nidaq card
disp('Initialising Nidaq card')
fs=400000;
ps.fs=fs;
ps.s = daq.createSession('ni');
ps.s.Rate = fs;
ao1=addAnalogOutputChannel(ps.s,'dev1',0,'Voltage'); % sound output
ao2=addAnalogOutputChannel(ps.s,'dev1',1,'Voltage'); % event output
ai1=addAnalogInputChannel(ps.s,'dev1',1,'Voltage'); % frame events input
ai2=addAnalogInputChannel(ps.s,'dev1',8,'Voltage'); % stimulus events input
ps.s.Channels(3).InputType = 'SingleEnded';
ps.s.Channels(4).InputType = 'SingleEnded';
ps.s.ExternalTriggerTimeout = 15;
ps.fs=ps.s.Rate; % check sample rate
if ps.fs~=fs
    disp('Sample rate not available!!!!')
    return
end

ps.fileIndex=1;
disp('Done');

%% prepare the stimulus

% Work out how long your recording will be - so we know how many frames to
% acquire
for ff = 1:length(ps.stimFiles)
    stim=[];
    stimInf = audioinfo([ps.stimFolder ps.stimFiles(ff).name]);
    if stimInf.SampleRate~=fs
        disp('STIM AT WRONG SAMPLE RATE!!')
        keyboard
    end
    ps.stimDur(ff) = stimInf.TotalSamples;
end

% pad with zeros at the end if not a whole number of seconds (easier if
% whole number)
wn = fs-mod(sum(ps.stimDur),fs);
if wn~=fs
    zeroPad = zeros(wn,2);
else
    zeroPad=[];
end
ps.nChunks = (sum(ps.stimDur)+length(zeroPad))/fs;

% Add trigger to start recording
triggerDuration = 500; % in samples
triggerAcquisition = [zeros(bkg*fs,1),[ones(triggerDuration,1)*5;zeros((bkg*fs)-triggerDuration,1)]]; % Initial trigger event to the 2P microscope

% Add the listeners for continuous playback/acquisition
ps.lh = addlistener(ps.s,'DataRequired',@presentStimCont_v3);
ps.la = addlistener(ps.s,'DataAvailable',@acquireCont_v2);
ps.fid = fopen(['E:\data\' ps.saveName],'a'); % open file for acquired data
ps.s.IsContinuous = true; % set nidaq to continuous mode
ps.counter=1; % start counters

stimD = sum(ps.stimDur)+length(zeroPad)+bkg*fs+5*fs; % total stimulus duration
disp(['stimulus duration: ' num2str(stimD/fs) ' s, ' num2str(stimD/fs/60) ' mins'])
fr = 29.874; % 1.25X %fr=29.543; % 1.5X
disp([num2str(fr) ' fps : ' num2str((stimD/fs)*fr) ' frames']); % estimate number of frames

%% present the stimulus
ps.firstChunk=1;
queueOutputData(ps.s,triggerAcquisition);
% Initialise the presentation/acquisition (the listeners take over after
% triggerAcquisition has been presented
ps.s.startBackground();
