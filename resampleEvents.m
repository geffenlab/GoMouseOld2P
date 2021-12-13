clear
stimLoc=['C:\Experiments\newData\'];
files = dir([stimLoc '*.txt']);
sForm = '%f%f';
nSamples = 4000000; % number of samples to read from text file in one go
fs=400000;
fid2 = fopen([stimLoc 'downsampled.txt'],'a');
for ii=1:length(files)
    fid = fopen([stimLoc files(ii).name]);
    frameDur=[]; eventsOn=[];eventsOff=[];
    tic
    while ~feof(fid)
        
        data = textscan(fid,sForm,nSamples,'Delimiter','\t');
        y = [data{1},data{2}];       
         y=resample(y,50000,400000);
        fprintf(fid2,'%f\t%f\n',y');
        toc
    end
end

fclose('all')