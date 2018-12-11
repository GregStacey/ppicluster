
% network noise / chromatogram noise
% dense/sparse

% TO DO:
%   1. In Craig netnoise, replace shufflenetwork with add+subtract.
%   1. In Craig, include all clustering algorithms.


%% 0. Get data

datadir = 'E:/Greg/PCP_SILAC/Runs/Craig/Final_interactome_13k/Output/tmp/';
datanames = {'Craig'};
tmp = dir([datadir 'data*rep*chan*.mat']);

fnsave = 'E:\Greg\ClusterReliable\data/data_save_08.mat';

java_path = 'E:/Greg/ClusterReliable/java/cluster_one-1.0.jar';

try load(fnsave)
catch
    
    %
    clear data
    data.file.name = cell(size(tmp));
    data.file.rep = nan(size(tmp));
    data.file.channel = nan(size(tmp));
    data.Chromatograms = cell(size(tmp));
    data.Protein = cell(size(tmp));
    for ii = 1:length(tmp)
        data.file.name{ii} = [tmp(ii).folder '\' tmp(ii).name];
        
        % meta info
        tmp2 = strsplit(tmp(ii).name, '_');
        data.file.rep(ii) = str2double(strrep(tmp2{2}, 'rep', ''));
        data.file.channel(ii) = str2double(strrep(strrep(tmp2{3}, 'chan', ''),'.mat',''));
        
        % data
        load(data.file.name{ii})
        data.Chromatograms{ii} = Chromatograms;
        data.Protein{ii} = Protein;
    end
    
    
    
    %% 1. Make chrom-noise score matrices
    
    % multiplicative noise, chrom_new = chrom * (1 + randn * s)
    %   s = [0 .1 .25 .5 1 2 5]
    
    sRange = [0 0.01 0.02 0.05 .1 0.15 .25 .5 1];
    data.chromnoise.Chromatograms = cell(length(data.file.name),length(sRange));
    data.chromnoise.score = cell(length(data.file.name),length(sRange));
    data.chromnoise.interactome = cell(length(data.file.name),length(sRange));
    for ii = 1:length(data.file.name)
        load(data.file.name{ii})
        for jj = 1:length(sRange)
            disp(['Making score matrix... ' num2str(ii) ' ' num2str(jj)])
            % add noise to chroms
            rnd = randn(size(data.Chromatograms{ii}));
            data.chromnoise.Chromatograms{ii,jj} = data.Chromatograms{ii}.*exp(rnd*sRange(jj));
            
            % pre-process chromatograms
            chromclean = data.chromnoise.Chromatograms{ii,jj};
            I = find(isnan(chromclean));
            chromclean(I)= 0.05 * rand(size(I));
            
            % calc new Dist matrices
            clear tmpdist
            tmpdist.Euc = squareform(pdist(chromclean,'euclidean'));              % 1. Euclidean distance
            tmpdist.R = -1*corr(chromclean');                                     % 2. Cleaned chromatogram R^2
            [R,p] = corrcoef(data.chromnoise.Chromatograms{ii,jj}','rows','pairwise');
            tmpdist.Rraw = 1 - R;                                                    % 3. Raw chromatogram R^2
            tmpdist.Rpraw = p;                                                       % 4. Raw chromatogram correlation p-value
            [~,mx] = max(chromclean,[],2);
            tmpdist.CoApex2 = squareform(pdist(mx,'euclidean'));                     % 6. Co-Apex score 2
            fn = fieldnames(tmpdist);
            NDistFields = length(fn);
            
            % calc score
            I = find(triu(ones(size(TP_Matrix)),1) == 1);
            possList = possibleInts(I(:));
            classList = TP_Matrix(I(:));
            fnames = fieldnames(tmpdist);
            DistList = nan(length(possList), length(fnames));
            for kk = 1:length(fnames)
                DistList(:,kk) = tmpdist.(fnames{kk})(I(:));
            end
            [scoreMatrix, feats_new] = scorenb(DistList,possList,classList);
            score = nanmedian(scoreMatrix,2);
            score(isnan(score)) = 0;
            
            % convert score to precision
            [x,I2] = sort(score, 'descend');
            score = score(I2);
            possList = possList(I2);
            classList = classList(I2);
            classList(possList==0) = nan;
            TP = cumsum(classList==1);
            FP = cumsum(classList==0);
            prec = TP ./ (TP + FP);
            prec(isnan(prec) | isinf(prec)) = 1;
            
            % turn score into a matrix
            [ia,ib] = ind2sub(size(TP_Matrix), I(I2));
            data.chromnoise.score{ii,jj} = nan(size(TP_Matrix));
            for kk = 1:length(ia)
                data.chromnoise.score{ii,jj}(ia(kk), ib(kk)) = prec(kk);
                data.chromnoise.score{ii,jj}(ib(kk), ia(kk)) = prec(kk);
            end
            
            % get 50% interactome
            [ia,ib] = find(data.chromnoise.score{ii,jj} >= 0.50);
            ibad = ia<ib;
            ia(ibad) = [];
            ib(ibad) = [];
            data.chromnoise.interactome{ii,jj} = [ia ib];
        end
    end
    
    
    
    %% 2a. Cluster data - chrom noise
    
    pars.best_p = 500;
    pars.best_dens = 0.1;
    pars.best_prec = 0.6;
    pars.best_I = 4;
    
    fnames1 = {'co_mcl' 'co' 'mcl'};
    fnames2 = {'mr' 'ga' 'sn' 'ppv' 'nmi' 'ari' 'coint' 'cocom'};
    for ii = 1:length(fnames1)
        data.chromnoise.(fnames1{ii}).cluster = cell(size(data.chromnoise.score));
        for jj = 1:length(fnames2)
            data.chromnoise.(fnames1{ii}).(fnames2{jj}) = nan(size(data.chromnoise.score));
        end
    end
    
    for ii = 1:length(data.file.name)
        % chrom_noise clusters
        for jj = 1:size(data.chromnoise.score,2)
            disp('Clustering net_noise')
            disp(['    chrom noise ' num2str(cRange(jj)*100) '%'])
            % make intMatrix
            intMatrix = data.chromnoise.score{ii,jj} > 0.5;
            intMatrix(intMatrix<0) = 0;
            intMatrix = (intMatrix - nanmin(intMatrix(:))) / (nanmax(intMatrix(:)) - nanmin(intMatrix(:)));
            if nansum(intMatrix(:))==0; continue; end
            
            for kk = 1:length(fnames1)
                % cluster
                if ismember('co_mcl',fnames1(kk))
                    disp('        co_mcl')
                    data.chromnoise.co_mcl.cluster{ii,jj} = clustone_mcl(intMatrix, pars, java_path);
                elseif ismember('co',fnames1(kk))
                    disp('        co')
                    data.chromnoise.co.cluster{ii,jj} = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path);
                elseif ismember('mcl',fnames1(kk))
                    disp('        mcl')
                    [~, data.chromnoise.mcl.cluster{ii,jj}] = mymcl(intMatrix, 1.8);
                end
                
                if isempty(data.chromnoise.(fnames1{kk}).cluster{ii,1}); continue; end
                if isempty(data.chromnoise.(fnames1{kk}).cluster{ii,jj}); continue; end
                
                % compare networks
                tmp = comparenetworks(data.chromnoise.(fnames1{kk}).cluster{ii,jj},...
                    data.chromnoise.(fnames1{kk}).cluster{ii,1},...
                    intMatrix,data.chromnoise.score{ii,1});
                fnames_compare = fieldnames(tmp);
                for mm = 1:length(fnames_compare)
                    data.chromnoise.(fnames1{kk}).(fnames_compare{mm})(ii,jj) = tmp.(fnames_compare{mm});
                end
                if ismember('mcl', fnames1{kk})
                    data.chromnoise.mcl.ari(ii,jj) = ari(data.chromnoise.mcl.cluster{ii,1},...
                        data.chromnoise.mcl.cluster{ii,jj});
                end
            end
        end
    end
    
    
    %% 2c. Cluster data - network noise
    
%     intMatrix = shufflenetwork(data.chromnoise.score{ii,1}, sRange(jj));
%     intMatrix = intMatrix > 0.5;
%     intMatrix(intMatrix<0) = 0;
%     intMatrix = (intMatrix - nanmin(intMatrix(:))) / (nanmax(intMatrix(:)) - nanmin(intMatrix(:)));
%     if nansum(intMatrix(:))==0; continue; end
    
    pars.best_p = 500;
    pars.best_dens = 0.1;
    pars.best_prec = 0.6;
    pars.best_I = 4;
    
    fnames0 = {'add' 'remove'};
    fnames1 = {'co_mcl' 'co' 'mcl'};
    fnames2 = {'mr' 'ga' 'sn' 'ppv' 'nmi' 'ari' 'coint' 'cocom'};
    for cc = 1:length(fnames0)
        for ii = 1:length(fnames1)
            data.netnoise.(fnames0{cc}).(fnames1{ii}).cluster = cell(size(data.chromnoise.score));
            for jj = 1:length(fnames2)
                data.netnoise.(fnames0{cc}).(fnames1{ii}).(fnames2{jj}) = nan(size(data.chromnoise.score));
            end
        end
    end
    
    for ii = 1:length(data.file.name)
        % chrom_noise clusters
        for jj = 1:size(data.chromnoise.score,2)
            disp('Clustering net_noise')
            for cc = 1:length(fnames0)
                % make intMatrix
                intMatrix = data.chromnoise.score{ii,1} > 0.5;
                if ismember(fnames0{cc}, 'add')
                    disp(['    add ' num2str(cRange(jj))])
                    intMatrix = addremovenetwork(intMatrix, 1 * cRange(jj));
                elseif ismember(fnames0{cc},'remove')
                    disp(['    remove ' num2str(cRange(jj))])
                    intMatrix = addremovenetwork(intMatrix, -1 * cRange(jj));
                end
                
                for kk = 1:length(fnames1)
                    % cluster
                    if ismember('co_mcl',fnames1(kk))
                        disp('        co_mcl')
                        data.netnoise.(fnames0{cc}).co_mcl.cluster{ii,jj} = clustone_mcl(intMatrix, pars, java_path);
                    elseif ismember('co',fnames1(kk))
                        disp('        co')
                        data.netnoise.(fnames0{cc}).co.cluster{ii,jj} = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path);
                    elseif ismember('mcl',fnames1(kk))
                        disp('        mcl')
                        [~, data.netnoise.(fnames0{cc}).mcl.cluster{ii,jj}] = mymcl(intMatrix, 1.8);
                    end
                    
                    
                    if isempty(data.netnoise.(fnames0{cc}).(fnames1{kk}).cluster{ii,1}); continue; end
                    if isempty(data.netnoise.(fnames0{cc}).(fnames1{kk}).cluster{ii,jj}); continue; end
                    
                    % compare networks
                    tmp = comparenetworks(data.netnoise.(fnames0{cc}).(fnames1{kk}).cluster{ii,jj},...
                        data.netnoise.(fnames0{cc}).(fnames1{kk}).cluster{ii,1},...
                        intMatrix,data.chromnoise.score{ii,1});
                    fnames_compare = fieldnames(tmp);
                    for mm = 1:length(fnames_compare)
                        data.netnoise.(fnames0{cc}).(fnames1{kk}).(fnames_compare{mm})(ii,jj) = tmp.(fnames_compare{mm});
                    end
                    if ismember('mcl', fnames1{kk})
                        data.netnoise.(fnames0{cc}).mcl.ari(ii,jj) = ari(data.netnoise.(fnames0{cc}).mcl.cluster{ii,1},...
                            data.netnoise.(fnames0{cc}).mcl.cluster{ii,jj});
                    end
                end
            end
        end
    end
    
    
    %% 3. Cluster corum
    
    pars.best_p = 500;
    pars.best_dens = 0.1;
    pars.best_prec = 0.6;
    pars.best_I = 4;
    
    % turn corum into a network
    fn = 'E:/Greg/ClusterReliable/data/allComplexes.txt';
    data.corum.network = corum2network(fn);
    
    cRange = sRange;
    for ii = 1:length(cRange)
        disp('Clustering corum... ')
        disp(['    add ' num2str(cRange(ii))])
        % ADD edges to corum-network
        intMatrix = addremovenetwork(data.corum.network, 1 * cRange(ii));
        
        disp('        clusterone+mcl')
        % cluster with clusterone+mcl
        data.corum.add.co_mcl.cluster{ii} = clustone_mcl(intMatrix, pars, java_path);
        % compare networks
        tmp = comparenetworks(data.corum.add.co_mcl.cluster{ii},data.corum.add.co_mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.add.co_mcl.mr(ii) = tmp.mr;
        data.corum.add.co_mcl.nmi(ii) = tmp.nmi;
        data.corum.add.co_mcl.ga(ii) = tmp.ga;
        data.corum.add.co_mcl.sn(ii) = tmp.sn;
        data.corum.add.co_mcl.ppv(ii) = tmp.ppv;
        data.corum.add.co_mcl.coint(ii) = tmp.coint;
        data.corum.add.co_mcl.cocom(ii) = tmp.cocom;
        
        disp('        clusterone')
        % cluster with clusterone
        data.corum.add.co.cluster{ii} = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path);
        % compare networks
        tmp = comparenetworks(data.corum.add.co.cluster{ii},data.corum.add.co.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.add.co.mr(ii) = tmp.mr;
        data.corum.add.co.nmi(ii) = tmp.nmi;
        data.corum.add.co.ga(ii) = tmp.ga;
        data.corum.add.co.sn(ii) = tmp.sn;
        data.corum.add.co.ppv(ii) = tmp.ppv;
        data.corum.add.co.coint(ii) = tmp.coint;
        data.corum.add.co.cocom(ii) = tmp.cocom;
        
        disp('        mcl')
        % cluster with mcl
        [~, data.corum.add.mcl.cluster{ii}] = mymcl(intMatrix, 1.8);
        % compare networks
        tmp = comparenetworks(data.corum.add.mcl.cluster{ii},data.corum.add.mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.add.mcl.mr(ii) = tmp.mr;
        data.corum.add.mcl.nmi(ii) = tmp.nmi;
        data.corum.add.mcl.ga(ii) = tmp.ga;
        data.corum.add.mcl.sn(ii) = tmp.sn;
        data.corum.add.mcl.ppv(ii) = tmp.ppv;
        data.corum.add.mcl.coint(ii) = tmp.coint;
        data.corum.add.mcl.cocom(ii) = tmp.cocom;
        data.corum.add.mcl.ari(ii) = ari(data.corum.add.mcl.cluster{ii},data.corum.add.mcl.cluster{1});
        
        
        % REMOVE edges to corum-network
        disp(['    remove ' num2str(cRange(ii))])
        intMatrix = addremovenetwork(data.corum.network, -1 * cRange(ii));
        
        disp('        clusterone+mcl')
        % cluster with clusterone+mcl
        data.corum.remove.co_mcl.cluster{ii} = clustone_mcl(intMatrix, pars, java_path);
        % compare networks
        tmp = comparenetworks(data.corum.remove.co_mcl.cluster{ii},data.corum.remove.co_mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.remove.co_mcl.mr(ii) = tmp.mr;
        data.corum.remove.co_mcl.nmi(ii) = tmp.nmi;
        data.corum.remove.co_mcl.ga(ii) = tmp.ga;
        data.corum.remove.co_mcl.sn(ii) = tmp.sn;
        data.corum.remove.co_mcl.ppv(ii) = tmp.ppv;
        data.corum.remove.co_mcl.coint(ii) = tmp.coint;
        data.corum.remove.co_mcl.cocom(ii) = tmp.cocom;
        
        disp('        clusterone')
        % cluster with clusterone
        data.corum.remove.co.cluster{ii} = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path);
        % compare networks
        tmp = comparenetworks(data.corum.remove.co.cluster{ii},data.corum.remove.co.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.remove.co.mr(ii) = tmp.mr;
        data.corum.remove.co.nmi(ii) = tmp.nmi;
        data.corum.remove.co.ga(ii) = tmp.ga;
        data.corum.remove.co.sn(ii) = tmp.sn;
        data.corum.remove.co.ppv(ii) = tmp.ppv;
        data.corum.remove.co.coint(ii) = tmp.coint;
        data.corum.remove.co.cocom(ii) = tmp.cocom;
        
        disp('        mcl')
        % cluster with mcl
        [~, data.corum.remove.mcl.cluster{ii}] = mymcl(intMatrix, 1.8);
        % compare networks
        tmp = comparenetworks(data.corum.remove.mcl.cluster{ii},data.corum.remove.mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.remove.mcl.mr(ii) = tmp.mr;
        data.corum.remove.mcl.nmi(ii) = tmp.nmi;
        data.corum.remove.mcl.ga(ii) = tmp.ga;
        data.corum.remove.mcl.sn(ii) = tmp.sn;
        data.corum.remove.mcl.ppv(ii) = tmp.ppv;
        data.corum.remove.mcl.coint(ii) = tmp.coint;
        data.corum.remove.mcl.cocom(ii) = tmp.cocom;
        data.corum.remove.mcl.ari(ii) = ari(data.corum.remove.mcl.cluster{ii},data.corum.remove.mcl.cluster{1});
    end
    
    save(fnsave, 'data');
end



