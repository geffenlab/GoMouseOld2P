% Finds all files recorded that day and extracts frame and event times - in
% samples, downsamples the data and re-saves
clear
stimLoc='C:\data\';
zipLoc = 'E:\data\';
% stimLoc=['C:\Users\labuser\Desktop'];

folders = dir(stimLoc);
folders(strcmp({folders.name},'.'))=[]; folders(strcmp({folders.name},'..'))=[];
today = datestr(now);
ind=1;
for ii=1:length(folders)
    dd = daysdif(folders(ii).date,today);
%     if dd<1
        f = dir([stimLoc folders(ii).name '\*.txt']);
        for jj=1:length(f)
            files{ind} = [stimLoc folders(ii).name '\' f(jj).name]; %#ok<*SAGROW>
            ind = ind+1;
        end
%     end
end



% Remove files you have already processed...
for ii=1:length(files)
    gr(ii) = any(strcmp([files{ii}(1:end-4) '_dsamp50kHz.txt'],files));
end
files(gr)=[];
files(contains(files,'_dsamp50kHz.txt'))=[];

% define some variables
sForm = '%f\t%f';
nSamples = 4000000; % number of samples to read from text file in one go
% fs=400000; % sample rate of recording
dfs = 50000; % sample rate to downsample to

% Check connected to data analysis pc
uigetdir('\\DESKTOP-GK8OVIP\data\','Check you are connected to the data analysis computer')

for ii=1:length(files)
    load([files{ii}(1:end-4) '_exptInfo.mat']);
    fs = exptInfo.fsStim;  % sample rate of recording/playback
    mouse = exptInfo.mouse;
    fid = fopen([files{ii}]);
    fid2 = fopen([files{ii}(1:end-4) '_dsamp50kHz.txt'],'a');
    frameOn=[]; eventOn=[];eventOff=[]; frameOff=[];
    tic
    ind=0;
    prevB = [];
    while ~feof(fid)
        
        data = textscan(fid,sForm,nSamples,'Delimiter','\t');
        y = [data{1},data{2}];
        if ~isempty(y)
            b=resample(y,50000,fs);
            fprintf(fid2,'%f\t%f\n',b');
        end
        
        y = [prevB;y]; 
        prevB = y(end-1:end,:);
        % Work out frame rate
        a = y(:,1)>1;
        dt = diff(a);
        ev = find(dt==1)+1;
        evOff = find(dt==-1)+1;
        
        frameOn = [frameOn;ev+(ind*nSamples)-1]; %#ok<*AGROW>
        frameOff = [frameOff;evOff+(ind*nSamples)-1];
        % Work out events
        a = y(:,2)>0.9;
        dt = diff(a);
        ev = find(dt==1)+1;
        evOff = find(dt==-1)+1;
        eventOn = [eventOn;ev+(ind*nSamples)-1];
        eventOff = [eventOff;evOff+(ind*nSamples)-1];
        toc
        ind=ind+1;
    end
    try
        save([files{ii}(1:end-4) '_exptInfo.mat'],'frameOn','frameOff','eventOff','eventOn','fs','-append')
    catch
        save([files{ii}(1:end-4) '_exptInfo.mat'],'frameOn','frameOff','eventOff','eventOn','fs')
    end
    fclose('all');
    zip([files{ii}(1:end-4) '_dsamp50kHz.zip'],[files{ii}(1:end-4) '_dsamp50kHz.txt'])
    movefile([files{ii}(1:end-4) '_dsamp50kHz.zip'],[zipLoc 'zippedEventTraces\'])
    us = strfind(files{ii},'_');
    aFold = ['\\DESKTOP-GK8OVIP\data\' mouse '\' datestr(now,'yyyymmdd') files{ii}(9:12) '_tifStacks'];
    if ~isdir(aFold); mkdir(aFold); end % make folder if does not exist
    copyfile([files{ii}(1:end-4) '_exptInfo.mat'],aFold);
%     if ~isdir([files{ii}(us(2)-4:us(2)-1) '\']); mkdir([files{ii}(us(2)-4:us(2)-1) '\']); end % make folder if does not exist
%     movefile([files{ii}(1:end-4) '_exptInfo.mat'],[files{ii}(us(2)-4:us(2)-1) '\']);
    delete([files{ii}(1:end-4) '_dsamp50kHz.txt']);
    delete([files{ii}]);
end

disp('finished extracting events');