function wavwrite_append(s, fn, chunk_size, fs, nbits)
%
%   s: stimulus
%   fn: file name
%   chunk_size:
%   fs:
%   nbits:

fid = [];
fmt = [];

if ~exist('chunk_size','var') || isempty(chunk_size)
    chunk_size = 1e6;
end

l = length(s);
n = floor(l / chunk_size);
last_chunk = l - (n * chunk_size); 

for i = 1:n + 1
    fprintf('\t\tWriting chunk %g/%g...\n',i,n+1);
    start_idx = (chunk_size * (i-1)) + 1;
    end_idx = start_idx + chunk_size - 1;
    switch(i)
        case 1
            write = 1;
        case n + 1
            write = 4;
            end_idx = l;
        otherwise
            write = 3;
    end
    [fid,fmt] = wavwriteStim(s(start_idx:end_idx,:),fs,nbits,fn,write,l,fid,fmt);
end