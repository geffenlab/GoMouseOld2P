function acquireCont(src, event)
global ps
if ps.counter < size(ps.stim,3)
%     disp([num2str(ps.counter) '/' num2str(size(ps.stim,3))])
%     ps.counter=ps.counter+1;
%     queueOutputData(ps.s,squeeze(ps.stim(:,:,ps.counter)));
disp('acquiring continuously')
a = ps.fs/10;
ps.data((ps.acqCounter-1)*a+1:ps.acqCounter*a,:) = event.Data;
ps.acqCounter=ps.acqCounter+1;
else
%     ps.s.stop();
%     delete(ps.la);
%     ps.s.IsContinuous = false;
end