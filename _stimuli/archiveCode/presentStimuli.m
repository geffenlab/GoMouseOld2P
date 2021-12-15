%% set variables
clear all
clear global
clc
delete(instrfindall)
global ps % this global is needed for the continuous nidaq presentation
repeats = 4;
bkg=5; % silence pre stim in s, make whole number
stimFile = '20170228_randReg_pool1585dB_1rep';
mouse = 'K055';
ps.saveName = [datestr(now,'yymmdd') mouse '_' stimFile '.mat'];

cd('C:\Users\Chris\Documents\GitHub\Kath\stimuli\')

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
stim=audioread([stimFile '.wav']);

stim=stim*10;
stim = repmat(stim,repeats,1); % number of repeats

% pad with zeros at the end if not a whole number of seconds
wn = round((1-mod(length(stim)/fs,1))*fs); 
if wn~=1
    stim = [stim;zeros(wn,2)];
end

% Add trigger to start recording
triggerDuration = 500; % in samples
stim = [[zeros(bkg*fs,1),[ones(triggerDuration,1)*5;zeros((bkg*fs)-triggerDuration,1)]];stim];

% Add some zeros at the end to make sure stimulus is not cut off
bkgPadding = zeros(fs*10,2);
stim = [stim;bkgPadding];

% Divide the stim into 1 second chunks for continuous playback using the
% nidaq
ps.stim = zeros(fs,2,length(stim)/fs);
for ii=1:length(stim)/fs
    ps.stim(:,:,ii) = stim((ii-1)*fs+1:ii*fs,:);
    disp([num2str(ii) '/' num2str(length(stim)/fs)])
end
clear stim bkgPadding

% Add the continuous playback and acquisition listeners
ps.lh = addlistener(ps.s,'DataRequired',@presentStimCont);
ps.la = addlistener(ps.s,'DataAvailable',@acquireCont);
ps.data=zeros(length(ps.stim(:))/2,2); % preallocate recording data
ps.s.IsContinuous = true;
ps.counter=1; % start counters
ps.acqCounter=1; % start counters

% display stimulus info
disp(['stimulus duration: ' num2str(length(ps.stim(:))/2/fs) ' s, ' num2str(length(ps.stim(:))/2/fs/60) ' mins'])
fr = 29.534;
disp([num2str(fr) ' fps : ' num2str(((length(ps.stim(:))/2)/fs)*fr) ' frames']); % estimate number of frames
%% present the stimulus
queueOutputData(ps.s,squeeze(ps.stim(:,:,1)));
ps.s.startBackground();  
% ps.data = [ps.data;ps.s.event.Data];
