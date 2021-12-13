% extract events from down sampled data...

% FIRST UNZIP YOUR FILE AND MOVE IT TO THE ZIPPED FILES FOLDER...

clear
stimLoc=['E:\data\zippedEventTraces\'];
files = dir([stimLoc '*.txt']);

% define some variables
sForm = '%f\t%f';
nSamples = 50000; % number of samples to read from text file in one go
fs=50000; % sample rate of recording

for ii=1:length(files)
    fid = fopen([stimLoc files(ii).name]);
    frameOn=[]; eventOn=[];eventOff=[]; frameOff=[];
    tic
    ind=0;
    prevB = [];
    z=[];
    counter = 0;
    while ~feof(fid)
        counter = counter+1;
        disp(counter)
        data = textscan(fid,sForm,nSamples,'Delimiter','\t');
        y = [data{1},data{2}];
        if ~isempty(y)
            b=resample(y,50000,400000);
        end
        
        y = [prevB;y];
        prevB = y(end-1:end,:);
        % Work out frame rate
        a = y(:,1)>1;
        dt = diff(a);
        ev = find(dt==1)+1;
        disp(['number of events ' num2str(length(ev))])
        nfr(counter) = length(ev);
        evOff = find(dt==-1)+1;

        frameOn = [frameOn;ev+(ind*nSamples)-1];
        frameOff = [frameOff;evOff+(ind*nSamples)-1];
        % Work out events
        a = y(:,2)>0.1;
        z = [z;a];
        dt = diff(a);
        ev = find(dt==1)+1;
        evOff = find(dt==-1)+1;
        eventOn = [eventOn;ev+(ind*nSamples)-1];
        eventOff = [eventOff;evOff+(ind*nSamples)-1];
        toc
        ind=ind+1;
    end
%     sf = input('enter expt file name: ');
%     save(sf,'frameOn','frameOff','eventOff','eventOn','fs','-append')
    fclose('all');
    pause()
end

disp('finished extracting events');