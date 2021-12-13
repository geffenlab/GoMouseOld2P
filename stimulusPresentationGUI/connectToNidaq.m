function s = connectToNidaq(fs,channelsIn,channelsOut)
% out speaker = ch 0
% out events = ch 1
% frame events IN = ch 1
% stimulus events IN = ch 8

% channelsIn = [0 8];
% channelsOut = [0 1 2];
% fs = 400000;

delete(instrfindall)
daqreset
% disp('Initialising NIDAQ card')
s = daq.createSession('ni');
s.Rate = fs;
for ii=1:length(channelsOut)
    addAnalogOutputChannel(s,'dev1',channelsOut(ii),'Voltage'); % sound output
end
for ii=1:length(channelsIn)
    addAnalogInputChannel(s,'dev1',channelsIn(ii),'Voltage'); % frame events input
    s.Channels(ii+length(channelsOut)).InputType = 'SingleEnded';
end
s.ExternalTriggerTimeout = 15;
s.TriggersPerRun = 1;
s.NotifyWhenScansQueuedBelow = 370000;
fs_check=s.Rate; % check sample rate
if fs~=fs_check
    disp('Sample rate not available!!!!')
    keyboard
    return
end

% disp('Connected to NIDAQ');
