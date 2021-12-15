function acquireCont_v2(src, event)
global ps
if mod(ps.counter,10)==0
    disp('acquiring continuously')
end
fprintf(ps.fid,'%f\t%f\n',event.Data');
