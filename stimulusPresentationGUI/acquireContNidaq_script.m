function acquireContNidaq_script(src, event)

global nc

    
chanIn = nc.nidaq.input;

if mod(nc.counter,10)==0
    disp('acquiring continuously')
end
if length(chanIn)==3
    fprintf(nc.fid,'%f\t%f\t%f\n',event.Data');
else
    fprintf(nc.fid,'%f\t%f\n',event.Data');
end