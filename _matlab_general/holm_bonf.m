function HB_sig = holm_bonf(p,alpha)
% actually this is the  Hochberg step-up procedure


m = length(p);
[p_sorted,p_index] = sort(p);


jj = alpha ./ ((m - (1:m))+1);

k = find(p_sorted<=jj,1,'last');
if isempty(k)
    k = 1;
end

H = [ones(1,k-1),zeros(1,m-(k-1))];

HB_sig = H(p_index);

