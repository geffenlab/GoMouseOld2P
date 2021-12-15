function presentStimCont_v3(hObject, eventdata)
%
global ps
if ps.firstChunk==1
    ps.ff=1; ps.jj=1; ps.rm=[]; ps.sv=0; 
    ps.firstChunk=0;
end
if ps.counter <= ps.nChunks
    disp(sprintf('%d/%d',ps.counter,ps.nChunks));
    nfc = floor((ps.stimDur(ps.ff)+ps.sv)/ps.fs);
    disp(ps.jj)
    if isempty(ps.rm)
        stim = audioread([ps.stimFolder ps.stimFiles(ps.ff).name],...
            [((ps.jj-1)*ps.fs+1)-ps.sv,ps.jj*ps.fs-ps.sv]); % read in 1 second chunks
        ps.jj=ps.jj+1;
    else
        indexing = [(ps.jj-1)*ps.fs+1,ps.jj*ps.fs-ps.sv];
        stim = audioread([ps.stimFolder ps.stimFiles(ps.ff).name],...
            indexing); % read in 1 second chunks
        ps.jj=ps.jj+1;
        ps.rm=[];
    end
    stim=[ps.rm;stim];
    stim = stim*10; % Get back to full level (.wav files are saved as stim/10 so need to *10)

    queueOutputData(ps.s,stim);
    ps.counter=ps.counter+1;
    
    if ps.jj>nfc
        x=mod((ps.stimDur(ps.ff)+ps.sv),ps.fs);
        if x~=0
            rm2 = audioread([ps.stimFolder ps.stimFiles(ps.ff).name],...
                [ps.stimDur(ps.ff)-x+1,ps.stimDur(ps.ff)]);
        else
            rm2=[];
        end
        ps.rm=rm2;
        ps.sv = length(ps.rm);
        ps.ff=ps.ff+1;
        ps.jj=1;
        if ps.ff>ps.nFiles
            endPadding = zeros(ps.fs*6,2);
            stim=[ps.rm*10;zeros(ps.fs-ps.sv,2);endPadding];
            queueOutputData(ps.s,stim);
            ps.counter=ps.counter+1;
        end
    end

else

    stop(ps.s);
    delete(ps.lh);
    delete(ps.la);
    fclose(ps.fid);
    fclose('all');
    % save everything
    exptInfo.mouse = ps.mouse;
    exptInfo.repeats = ps.reps;
    exptInfo.stimFiles = {ps.stimFiles.name};
    b = unique(exptInfo.stimFiles);
    for ii=1:length(b)
        a = load([ps.stimFolder b{ii}(1:end-4) '_stimInfo.mat']);   
        exptInfo.stimInfo{ii} = a.stimInfo;
    end
    exptInfo.preStimSilence = ps.preStimSil;
    exptInfo.fsStim = ps.fs;
    exptInfo.yaw = ps.yaw;
    exptInfo.pitch = ps.pitch;
    save(['E:\data\' ps.saveName(1:end-4) '_exptInfo.mat'],'exptInfo')
    disp('FINISHED PRESENTING')
end



