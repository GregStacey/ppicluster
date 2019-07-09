function avgsilscore = mysilhouette(clusts, intMatrix)

% convert intMatrix into a distance matrix
intMatrix(eye(size(intMatrix))==1) = nan;
intMatrix = 1 - intMatrix;

% turn clusts (cell array) into 2d list of cluster ids
clustID = nan(size(intMatrix,1), 100);
%clustID = nan(size(intMatrix,1), 1);
for ii = 1:length(clusts)
    for jj = 1:length(clusts{ii})
        %clustID{clusts{ii}(jj)} = [clustID{clusts{ii}(jj)} ii];
        this_row = clusts{ii}(jj);
        nn = min(find(isnan(clustID(this_row,:))));
        clustID(this_row, nn) = ii;
    end
end
unqclusts = 1:length(clusts);

% give each data point a silhouette score
silscore = nan(size(intMatrix,1),1);
for ii = 1:length(silscore)
    all_clusts = clustID(ii, find(not(isnan(clustID(ii,:)))));
    if isempty(all_clusts); continue; end
    this_clust = randsample(all_clusts, 1);

    ia = find(clustID == this_clust);
    [ia, ~] = ind2sub(size(clustID), ia);
    ai = nanmean(nanmean(intMatrix(ia,ii)));
    
    d_other_clusts = nan(length(unqclusts), 1);
    for jj = 1:length(d_other_clusts)
        that_clust = unqclusts(jj);
        ib = find(clustID == that_clust);
        [ib, ~] = ind2sub(size(clustID), ib);
        d_other_clusts(jj) = nanmean(nanmean(intMatrix(ib,ii)));
    end
    d_other_clusts(all_clusts) = nan;
    bi = nanmin(d_other_clusts);
    
    silscore(ii) = (bi - ai) / max(ai, bi);
end

avgsilscore = nanmean(silscore);

end