function presentStimCont(hObject, eventdata, handles)
global ps
if ps.counter < size(ps.stim,3)
    ps.counter=ps.counter+1;
    disp([num2str(ps.counter) '/' num2str(size(ps.stim,3))])
    queueOutputData(ps.s,squeeze(ps.stim(:,:,ps.counter)));
else
    ps.fileIndex=ps.fileIndex+1;
    
%    pause(1.1);
     ps.s.stop();
    delete(ps.lh);
    delete(ps.la);
    fclose(ps.fid);
    ps.s.IsContinuous = false;
    pause(1.1)
%     y = ps.data;
%     save(['C:\Experiments\newData\' ps.saveName],'y','-v7.3');
end
