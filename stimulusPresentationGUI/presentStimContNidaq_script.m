function presentStimContNidaq_script(src, event)
%
disp('Looking for new stim')
global nc
% try
%     [~,chanOut] = getNidaqSettings(handles);
% catch
    chanOut = nc.nidaq.output;
% end

if nc.firstChunk==1
    nc.ff = 1; nc.jj = 1; nc.rm = []; nc.sv = 0;
    nc.firstChunk = 0; 
end

if nc.counter <= nc.nChunks
%     sprintf('%d/%d\n',nc.counter,nc.nChunks);
    disp('Loading next bit of stimulus')
    nfc = floor((nc.stimDur(nc.ff)+nc.sv)/nc.fs);
    disp([num2str(nc.jj) '/' num2str(nfc)])
    if isempty(nc.rm)
        stim = audioread(nc.stimFiles{nc.ff},...
            [((nc.jj-1)*nc.fs+1)-nc.sv,nc.jj*nc.fs-nc.sv]); % read in 1 second chunks
        if length(chanOut)==3
            stim(:,3) = audioread('temp_motionCamTrigger.wav',...
                [((nc.jj-1)*nc.fs+1)-nc.sv,nc.jj*nc.fs-nc.sv]); % read in 1 second chunks
        end
        nc.jj=nc.jj+1;
        disp('Stim loaded - 1')
    else
        indexing = [(nc.jj-1)*nc.fs+1,nc.jj*nc.fs-nc.sv];
        stim = audioread(nc.stimFiles{nc.ff},indexing); % read in 1 second chunks
        if length(chanOut)==3
            stim(:,3) = audioread('temp_motionCamTrigger.wav',indexing);
        end
        nc.jj = nc.jj+1;
        nc.rm = [];
         disp('Stim loaded - 2')
    end
    stim = [nc.rm;stim];
    stim = stim*10; % Get back to full level (.wav files are saved as stim/10 so need to *10)
    
    queueOutputData(nc.s,stim);
    nc.counter = nc.counter+1;
    disp('Queued output data')
    disp(['Nidaq running? ' num2str(nc.s.IsRunning)])
    
    if nc.jj>nfc 
        disp('Checking if any stim remaining')
        x=mod((nc.stimDur(nc.ff)+nc.sv),nc.fs);
        if x~=0
            rm2 = audioread(nc.stimFiles{nc.ff},...
                [nc.stimDur(nc.ff)-x+1,nc.stimDur(nc.ff)]);
        else
            rm2=[];
        end
        nc.rm=rm2;
        nc.sv = length(nc.rm);
        nc.ff = nc.ff+1;
        nc.jj=1;
        if nc.ff > nc.nFiles
            disp('Finishing up with a few zeros... wait for it!')
            endPadding = zeros(nc.fs*6,2);
            stim=[nc.rm*10;zeros(nc.fs-nc.sv,2);endPadding];
            if length(chanOut)==3 % add in the motion cammera
                pulse = [ones(0.001*nc.fs,1)*3;zeros(0.049*nc.fs,1)]; % 20 Hz frame rate
                dur = round(length(stim)/nc.fs);
                stim(:,3) = repmat(pulse,20*dur,1);
            end            
            queueOutputData(nc.s,stim);
            nc.counter=nc.counter+1;
        end
    end
    
else
    queueOutputData(nc.s,[0,0]);
    disp('End of stim')
    nc.blockN = nc.blockN+1;
    stop(nc.s);
    delete(nc.lh);
    if isfield(nc,'la')
        delete(nc.la);
    end
    fclose(nc.fid);
    fclose('all');
    % save everything
%     set(handles.text35,'String',['Saving block ' num2str(nc.blockN)])
    exptInfo.mouse = nc.mouse;
    exptInfo.stimFiles = nc.stimFiles;
    b = unique(exptInfo.stimFiles);
    for ii=1:length(b)
        a = load([b{ii}(1:end-4) '_stimInfo.mat']);
        exptInfo.stimInfo{ii} = a;%.stimInfo;
    end
    exptInfo.preStimSilence = nc.preStimSil;
    exptInfo.fsStim = nc.fs;
    exptInfo.yaw = nc.yaw;
    exptInfo.pitch = nc.pitch; %#ok<STRNU>
    save([nc.fn(1:end-4) '_exptInfo.mat'],'exptInfo')
    disp('stuff saved - FINISHED PRESENTING')
%      set(handles.text35,'String',['Block ' num2str(nc.blockN-1) ' of ' num2str(nc.nBlocks) ' saved'])
    
    %% NOW ADD IN A NEW FUNCTION, 'PLAYNEXTBLOCK'
%     if nc.blockN<=nc.nBlocks
%         playNextBlock(handles);
%     end
%     
end



