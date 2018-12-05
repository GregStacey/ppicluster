function Members2 = clustone_mcl(intMatrix, pars)

% Cluster intmatrix with ClusterONE+MCL

Members2 = [];
best_p = pars.best_p;
best_dens = pars.best_dens;
best_prec = pars.best_prec;
best_I = pars.best_I;

% ClusterONE - round 1
COne_members = myclusterone(intMatrix, best_p, 0);
if isempty(COne_members); return; end

% MCL - round 1
inflate_param = best_I;
cc = 0;
MCL_members = cell(10^5,1);
for kk = 1:length(COne_members)
    tmp = intMatrix(COne_members{kk},COne_members{kk});
    tmp2 = mymcl(tmp, inflate_param);
    Icomps = find(sum(tmp2,2)>2);
    for uu = 1:length(Icomps)
        cc = cc+1;
        MCL_members{cc} = COne_members{kk}(tmp2(Icomps(uu),:)>0);
    end
end
MCL_members = MCL_members(1:cc);

Members2 = MCL_members;
dens_comp = nan(size(Members2));
for kk = 1:size(Members2)
    I = Members2{kk};
    m = intMatrix(I,I);
    n = length(I);
    dens_comp(kk) = sum(m(:)) / (n * (n-1)/2);
    if n<3 || dens_comp(kk) < best_dens
        Members2{kk} = [];
    end
end
Members2(cellfun('isempty',Members2)) = [];
if isempty(Members2); return; end

% Run the remaining interactions through MCL
Iincomplex = unique([Members2{:}]);
intMatrix2 = intMatrix;
intMatrix2(Iincomplex,Iincomplex)=0;
% need to remove columns/rows with 0 interactions
Inot0 = find(nansum(intMatrix2)>0);
intMatrix2 = intMatrix2(Inot0,Inot0);
tmp2 = mymcl(intMatrix2, inflate_param);
Icomps = find(sum(tmp2,2)>2);
cc = length(Members2);
for uu = 1:length(Icomps)
    cc = cc+1;
    Members2{cc} = Inot0(tmp2(Icomps(uu),:)>0);
end

% Final removal of low density complexes
dens_comp = nan(size(Members2));
for kk = 1:size(Members2,1)
    I = Members2{kk};
    m = intMatrix(I,I);
    n = length(I);
    dens_comp(kk) = sum(m(:)) / (n * (n-1));
    if n<3 || dens_comp(kk) < best_dens
        Members2{kk} = [];
    end
end
Members2(cellfun('isempty',Members2)) = [];
if isempty(Members2); return; end

% Merge complexes with Jaccard>.9
if not(isequal(Members2{1}, Members2))
    JJ = nan(length(Members2), length(Members2));
    for kk = 1:length(Members2)
        for rr = 1:length(Members2)
            if kk>=rr; continue; end
            JJ(rr,kk) = length(intersect(Members2{kk}, Members2{rr})) / ...
                length(unique( [Members2{kk},Members2{rr}] ));
            if JJ(rr,kk)==1
                Members2{kk} = [];
                JJ(rr,kk) = 0;
            end
        end
    end
    % Remove empty complexes
    I = cellfun('isempty',Members2);
    Members2(I) = [];
end