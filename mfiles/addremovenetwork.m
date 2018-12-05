function M1 = addremovenetwork(M0,ff)
% noise network by adding OR removing 
%   M0 = un-noised network, square connection matrix
%   ff = fraction of indices to add OR remove. ff>0 means addition, ff<0
%        means removal.
%
% NB this ONLY works for SPARSE networks.

if sum(not(M0(:)==0)) + sum(not(M0(:)==1)) < length(M0(:))
    error('Only sparse binary networks with values of 0 and 1, please.')
end

if ff==0
    M1 = M0; 
    return 
end

% get un-noised network
I0 = find(triu(M0)==1);
mm = size(M0,1);

% how many edges to add/remove?
nn = round(ff * length(I0));

% make network-noise matrix (triu)
I0_noise = [];
if nn<0
    % remove edges
    I0_noise = randsample(I0, length(I0) + nn);
elseif nn>0
    % add edges
    Itriu = find(triu(ones(mm,mm),1)==1 & not(M0==1));
    Iadd = randsample(Itriu, nn);
    I0_noise = [I0; Iadd];
end
M1 = zeros(mm,mm);
M1(I0_noise) = 1;

% make matrix symmetric
M1 = tril(M1') + triu(M1);
