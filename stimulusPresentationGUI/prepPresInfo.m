function [presInfo] = prepPresInfo(handles)
global pm nc
% Get info to save to recording file
fs = str2double(get(handles.edit9,'String'));
presInfo.fs = fs;
presInfo.FILT = get(handles.edit10,'String');
presInfo.recType = get(handles.text17,'String');
oz = get(handles.text34,'String');
presInfo.opticalZoom = str2double(oz(1:end-1));
presInfo.mouse = pm.mouse;
presInfo.pitch = str2double(get(handles.edit2,'String'));
presInfo.yaw = str2double(get(handles.edit1,'String'));
pitch = presInfo.pitch; yaw = presInfo.yaw;
if pitch~=0 || yaw~=0
    save([pm.mouseFolder pm.mouse '\' pm.mouse 'pitchYaw.mat'],'pitch','yaw');
end

nc.pitch = presInfo.pitch;
nc.yaw = presInfo.yaw;

% Prepare info for presenting the stimuli
contents = cellstr(get(handles.listbox3,'String'));
sf = strcat(pm.wavFolders,contents);
d = get(handles.uitable2,'Data');
if ~iscell(d)
    d = num2cell(d); % blocks and repeats
end
if ~isempty(d)
    blocks = [d{:,1}];
    reps = [d{:,2}];
    stimFiles = cell(sum(reps),1); 
    b = zeros(sum(reps),1);
    ind = 1;
    if any(reps>1)
        for ii=1:length(reps)
            for jj=1:reps(ii)
                stimFiles(ind,1) = sf(ii);
                b(ind,1) = blocks(ii);
                ind = ind+1;
            end
        end
    else
        stimFiles = sf;
        b = blocks;
    end
    
    % sort by block
    [b,blockOrder] = sort(b);
    stimFiles = stimFiles(blockOrder);
    presInfo.blocks = b;
    presInfo.nBlocks = length(unique(b));
    presInfo.preStimSil = str2double(get(handles.edit4,'String'));
    presInfo.stimFiles = stimFiles;
    
    % Work out how long your recording will be - so we know how many frames to
    % acquire
    stimDur = cell(1,presInfo.nBlocks);
    for bb=1:presInfo.nBlocks
        ind = find(b==bb);
        for ff = 1:length(ind)
            stimInf = audioinfo(stimFiles{ind(ff)});
            set(handles.edit9,'String',num2str(stimInf.SampleRate))
            fs = str2double(get(handles.edit9,'String'));
            presInfo.fs = fs;
            if stimInf.SampleRate~=fs
                disp('STIM AT WRONG SAMPLE RATE!!')
                keyboard
            end
            stimDur{bb}(ff) = stimInf.TotalSamples;
        end
        
        
        %  pad with zeros at the end if not a whole number of seconds for each block
        % (easier if whole number)
        wn = fs-mod(sum(stimDur{bb}),fs);
        if wn~=fs
            zeroPad = zeros(wn,2);
        else
            zeroPad=[];
        end
        
        presInfo.nChunks(bb) = (sum(stimDur{bb})+length(zeroPad))/fs;
        
    end
    
    presInfo.stimDur = stimDur;
    
    [~,chanOut] = getNidaqSettings(handles);
    
    % Add trigger to start recording
    triggerDuration = 0.1*fs; % in samples
    presInfo.triggerAcquisition = [zeros(presInfo.preStimSil*fs,1),...
        [ones(triggerDuration,1)*5;zeros((presInfo.preStimSil*fs)...
        -triggerDuration,1)]]; % Initial trigger event to the 2P microscope
    
    if length(chanOut)==3 && chanOut(3)==2
        pulse = [ones(0.001*fs,1)*3;zeros(0.049*fs,1)]; % 20 Hz frame rate
        dur = round(length(presInfo.triggerAcquisition)/fs);
        presInfo.triggerAcquisition(:,3) = repmat(pulse,20*dur,1);
    elseif length(chanOut)==3 && chanOut(3)==3
        presInfo.triggerAcquisition(:,3) = zeros(length(presInfo.triggerAcquisition),1);
    elseif length(chanOut)==4
        presInfo.triggerAcquisition(:,3:4) = zeros(length(presInfo.triggerAcquisition),2);
    end
        
    
    presInfo.stimD = ceil((sum(presInfo.nChunks)*fs+...
        (presInfo.preStimSil*fs*presInfo.nBlocks)+...
        5*fs*presInfo.nBlocks)/fs); % total stimulus duration
    set(handles.text18,'String',sprintf('%0.1f',presInfo.stimD/60));
    
    
%     if length(chanOut)==3
% %         pulse = [ones(0.001*fs,1)*3;zeros(0.049*fs,1)]; % 20 Hz frame rate
%         dur = round(presInfo.nChunks+2);
%         motionCamTriggers = repmat(pulse,20*dur,1)/10;
%         audiowrite('temp_motionCamTrigger.wav',motionCamTriggers,fs,...
%             'bitsPerSample',16);
%     end
    
    
    
%     if strcmp(get(handles.text17,'String'),'2P imaging')
        fr = str2double(get(handles.edit6,'String'));
        set(handles.text19,'String',num2str(ceil(presInfo.stimD*fr)));
        %     disp([num2str(fr) ' fps : ' num2str((stimD/fs)*fr) ' frames']); % estimate number of frames
%     end
    
end