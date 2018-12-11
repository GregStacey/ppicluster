function output = comparenetworks(network, ref_network, matrix, ref_matrix)

% matching ratio
output.mr = matchingratio(ref_network, network);

% geometric accuracy
[output.ga,output.sn,output.ppv] = geomacc(ref_network, network);

% normalized mutual information
output.nmi = nmi(ref_network, network);

% co-complex probability
iterMax = 1000;
tmp = nan(iterMax, 1);
kk = 0;
while kk<=iterMax
    % pick a starting pair
    ia = randsample(length(ref_network), 1);
    ib1 = randsample(length(ref_network{ia}), 1);
    ib2 = randsample(length(ref_network{ia}), 1);
    if ib1==ib2; continue; end
    prot1 = ref_network{ia}(ib1);
    prot2 = ref_network{ia}(ib2);
    
    % are the in the target complex set?
    ia1 = cellfun(@(x) ismember(x, prot1), network, 'UniformOutput', 0);
    ia2 = cellfun(@(x) ismember(x, prot2), network, 'UniformOutput', 0);
    ia1 = find(cellfun(@(x) nansum(x)>0, ia1));
    ia2 = find(cellfun(@(x) nansum(x)>0, ia2));
    kk = kk+1;
    tmp(kk) = ~isempty(intersect(ia1, ia2));
end
output.cocom = nanmean(tmp);


% co-interactome probability
output.coint = nan;
try
    % reference 50% interactome
    [ia,ib] = find(ref_matrix >= 0.50);
    ibad = ia<ib;
    ia(ibad) = [];
    ib(ibad) = [];
    ref_interactome = [ia ib];
    
    % network 50% interactome
    [ia,ib] = find(matrix >= 0.50);
    ibad = ia<ib;
    ia(ibad) = [];
    ib(ibad) = [];
    interactome = [ia ib];
    
    iterMax = 1000;
    tmp = zeros(iterMax, 1);
    kk = 0;
    while kk<=iterMax
        % pick a starting pair
        ia = randsample(length(ref_interactome), 1);
        prot1 = ref_interactome(ia,1);
        prot2 = ref_interactome(ia,2);
        
        % are they in the target interactome?
        kk = kk+1;
        tmp(kk) = sum((ismember(interactome(:,1), prot1) & ismember(interactome(:,2), prot2)) ...
            | (ismember(interactome(:,2), prot1) & ismember(interactome(:,1), prot2))) > 0;
    end
    output.coint = nanmean(tmp);
end
