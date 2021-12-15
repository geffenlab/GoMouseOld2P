function presentStimCont_v2(hObject, eventdata)
global ps
if ps.counter < size(ps.stim,3)
    ps.counter=ps.counter+1;
    disp([num2str(ps.counter) '/' num2str(size(ps.stim,3))])
    queueOutputData(ps.s,squeeze(ps.stim(:,:,ps.counter)));
else
    ps.fileIndex=ps.fileIndex+1;
    
    if ps.fileIndex<=ps.nFiles % If still more files to play
        
    else % Stop the playback if played all files
        pause(1.1);
        ps.s.stop();
        delete(ps.lh);
        delete(ps.la);
        fclose(ps.fid);
        ps.s.IsContinuous = false;
    end

end
