function M1 = shufflenetwork(M0,ff)
% noise network by shuffling indices of triu(M0)
%   M0 = un-noised network, square connection matrix
%   ff = fraction of indices to shuffle
%
% NB this is designed for a DENSE network, and it will work only on average 
% with sparse binary networks.

% get un-noised network
mm = size(M0,1);

% how many edges to shuffle?
nn = round(ff * mm*(mm-1)/2);

% find nn indices in the upper right
Itriu = find(triu(ones(mm,mm),1)==1);
Ishuffle = randsample(Itriu, nn);

% make network-noise matrix (triu)
if length(Ishuffle)>2
    Iinsert = [Ishuffle(end); Ishuffle(1:end-1)];
else
    Iinsert = [];
end
M1 = M0;
M1(Ishuffle) = M0(Iinsert);

% make matrix symmetric
M1 = tril(M1') + triu(M1);
