
% network noise / chromatogram noise
% dense/sparse

% TO DO:
%   1. In Craig netnoise, replace shufflenetwork with add+subtract.
%   1. In Craig, include all clustering algorithms.


%% 0. Get data

clear datadir
datadir{1} = 'E:/Greg/PCP_SILAC/Runs/Craig/Final_interactome_13k/Output/tmp/';
datadir{2} = 'E:/Greg/PCP_SILAC/Runs/Nick_Apoptosis/Output/tmp/';
datadir{3} = 'E:/Greg/PCP_SILAC/Runs/Anders_2012/Prince_paper_run/Output/tmp/';
datadir{4} = 'E:/Greg/PCP_SILAC/Runs/Nick_HeLa/Output/tmp/';

cc = 0;
for ii = 1:length(datadir)
    x = dir([datadir{ii} 'data*rep*chan*.mat']);
    for jj = 1:length(x)
        cc = cc+1;
        data.file.name{cc} = [x(jj).folder '/' x(jj).name];
    end
end

fnsave = 'E:\Greg\ClusterReliable\data/data_save_10.mat';

java_path = 'E:/Greg/ClusterReliable/java/cluster_one-1.0.jar';

try load(fnsave)
catch
    
    %
    clear data
    data.file.name = cell(size(data.file.name));
    data.file.rep = nan(size(data.file.name));
    data.file.channel = nan(size(data.file.name));
    data.Chromatograms = cell(size(data.file.name));
    data.Protein = cell(size(data.file.name));
    for ii = 21:length(data.file.name)
        % meta info
        fn = strsplit(data.file.name{ii}, '\');
        fn = fn{end};
        tmp2 = strsplit(fn, '_');
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
    
    data.chromnoiseRange = [0 0.01 0.02 0.05 .1 0.15 .25 .5 1];
    data.chromnoise.Chromatograms = cell(length(data.file.name),length(data.chromnoiseRange));
    data.chromnoise.score = cell(length(data.file.name),length(data.chromnoiseRange));
    data.chromnoise.interactome = cell(length(data.file.name),length(data.chromnoiseRange));
    for ii = 21:length(data.file.name)
        load(data.file.name{ii})
        for jj = 1:length(data.chromnoiseRange)
            disp(['Making score matrix... ' num2str(ii) ' ' num2str(jj)])
            % add noise to chroms
            rnd = randn(size(data.Chromatograms{ii}));
            data.chromnoise.Chromatograms{ii,jj} = data.Chromatograms{ii}.*exp(rnd*data.chromnoiseRange(jj));
            
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
            disp(['N predicted interactions: ' num2str(length(data.chromnoise.interactome{ii,jj}))])
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
            disp('Clustering chrom_noise')
            disp(['    chrom noise ' num2str(data.chromnoiseRange(jj)*100) '%'])
            % make intMatrix
            intMatrix = data.chromnoise.score{ii,jj} > 0.7;
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
                if ismember('mcl', fnames1(kk))
                    data.chromnoise.mcl.ari(ii,jj) = ari(data.chromnoise.mcl.cluster{ii,1},...
                        data.chromnoise.mcl.cluster{ii,jj});
                end
            end
        end
    end
    
    
    %% 2c. Cluster data - network noise
    
    pars.best_p = 500;
    pars.best_dens = 0.1;
    pars.best_prec = 0.6;
    pars.best_I = 4;
    
    fnames0 = {'add' 'remove' 'shuffle'};
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
    
    data.addnoiseRange = [0 0.01 0.05 .25 .5 1 2 3 5];
    data.remnoiseRange = [0 0.01 0.02 0.05 .1 0.15 .25 .5 0.95];
    data.shufflenoiseRange = [0 0.01 0.02 0.05 .1 0.15 .25 .5 1];
    for ii = 1:length(data.file.name)
        % chrom_noise clusters
        for jj = 1:size(data.chromnoise.score,2)
            disp('Clustering net_noise')
            for cc = 1:length(fnames0)
                % make intMatrix
                intMatrix = data.chromnoise.score{ii,1} > 0.5;
                if ismember(fnames0{cc}, 'add')
                    disp(['    add ' num2str(data.addnoiseRange(jj))])
                    intMatrix = addremovenetwork(intMatrix, 1 *  data.addnoiseRange(jj));
                elseif ismember(fnames0{cc},'remove')
                    disp(['    remove ' num2str(data.remnoiseRange(jj))])
                    intMatrix = addremovenetwork(intMatrix, -1 * data.remnoiseRange(jj));
                elseif ismember(fnames0{cc},'shuffle')
                    disp(['    shuffle ' num2str(data.remnoiseRange(jj))])
                    intMatrix = shufflenetwork(intMatrix, data.shufflenoiseRange(jj));
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
                    if ismember('mcl', fnames1(kk))
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
    [data.corum.network, data.corum.proteins] = corum2network(fn);
    
    for ii = 1:length(data.addnoiseRange)
        disp('Clustering corum... ')
        disp(['    add ' num2str(data.addnoiseRange(ii))])
        % ADD edges to corum-network
        intMatrix = addremovenetwork(data.corum.network, 1 * data.addnoiseRange(ii));
        
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
        disp(['    remove ' num2str(data.remnoiseRange(ii))])
        intMatrix = addremovenetwork(data.corum.network, -1 * data.remnoiseRange(ii));
        
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
        
        
        
        % SHUFFLE edges of corum-network
        disp(['    shuffle ' num2str(data.shufflenoiseRange(ii))])
        intMatrix = shufflenetwork(data.corum.network, data.shufflenoiseRange(ii));
        
        disp('        clusterone+mcl')
        % cluster with clusterone+mcl
        data.corum.shuffle.co_mcl.cluster{ii} = clustone_mcl(intMatrix, pars, java_path);
        % compare networks
        tmp = comparenetworks(data.corum.shuffle.co_mcl.cluster{ii},data.corum.shuffle.co_mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.shuffle.co_mcl.mr(ii) = tmp.mr;
        data.corum.shuffle.co_mcl.nmi(ii) = tmp.nmi;
        data.corum.shuffle.co_mcl.ga(ii) = tmp.ga;
        data.corum.shuffle.co_mcl.sn(ii) = tmp.sn;
        data.corum.shuffle.co_mcl.ppv(ii) = tmp.ppv;
        data.corum.shuffle.co_mcl.coint(ii) = tmp.coint;
        data.corum.shuffle.co_mcl.cocom(ii) = tmp.cocom;
        
        disp('        clusterone')
        % cluster with clusterone
        data.corum.shuffle.co.cluster{ii} = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path);
        % compare networks
        tmp = comparenetworks(data.corum.shuffle.co.cluster{ii},data.corum.shuffle.co.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.shuffle.co.mr(ii) = tmp.mr;
        data.corum.shuffle.co.nmi(ii) = tmp.nmi;
        data.corum.shuffle.co.ga(ii) = tmp.ga;
        data.corum.shuffle.co.sn(ii) = tmp.sn;
        data.corum.shuffle.co.ppv(ii) = tmp.ppv;
        data.corum.shuffle.co.coint(ii) = tmp.coint;
        data.corum.shuffle.co.cocom(ii) = tmp.cocom;
        
        disp('        mcl')
        % cluster with mcl
        [~, data.corum.shuffle.mcl.cluster{ii}] = mymcl(intMatrix, 1.8);
        % compare networks
        tmp = comparenetworks(data.corum.shuffle.mcl.cluster{ii},data.corum.shuffle.mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.shuffle.mcl.mr(ii) = tmp.mr;
        data.corum.shuffle.mcl.nmi(ii) = tmp.nmi;
        data.corum.shuffle.mcl.ga(ii) = tmp.ga;
        data.corum.shuffle.mcl.sn(ii) = tmp.sn;
        data.corum.shuffle.mcl.ppv(ii) = tmp.ppv;
        data.corum.shuffle.mcl.coint(ii) = tmp.coint;
        data.corum.shuffle.mcl.cocom(ii) = tmp.cocom;
        data.corum.shuffle.mcl.ari(ii) = ari(data.corum.shuffle.mcl.cluster{ii},data.corum.shuffle.mcl.cluster{1});
    end
    
    
    %% 4. Cluster corum with "bad" parameters
    
    pars.best_p = 2;
    pars.best_dens = 0.3;
    pars.best_prec = 0.6;
    pars.best_I = 2;
    
    for ii = 1:length(data.addnoiseRange)
        % SHUFFLE edges of corum-network
        disp(['    shuffle ' num2str(data.shufflenoiseRange(ii))])
        intMatrix = shufflenetwork(data.corum.network, data.shufflenoiseRange(ii));
        
        disp('        clusterone+mcl')
        % cluster with clusterone+mcl
        data.corum.shuffle_bad.co_mcl.cluster{ii} = clustone_mcl(intMatrix, pars, java_path, 1);
        
        disp('        clusterone')
        % cluster with clusterone
        data.corum.shuffle_bad.co.cluster{ii} = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path, 1);
        
        disp('        mcl')
        % cluster with mcl
        [~, data.corum.shuffle_bad.mcl.cluster{ii}] = mymcl(intMatrix, 4);
    end
    
    
    %% 5. choose parameters to make large clusters
    intMatrix = data.corum.network;
    
    % co
    pRange = [1, 50, 100, 500, 5000];
    densRange = [0, 0.1, 0.2, 0.3, 0.4];
    size_co = nan(length(pRange), length(densRange));
    for ii = 1:length(pRange)
        for jj = 1:length(densRange)
            clusts = clusterone_java(intMatrix, pRange(ii), densRange(jj), java_path);
            nn = nan(size(clusts));
            for uu = 1:length(nn)
                nn(uu) = length(clusts{uu});
            end
            nn = nn(nn>=3);
            size_co(ii,jj) = median(nn);
        end
    end

    
    % mcl
    Irange = [2, 4, 8, 16, 32];
    Erange = [2, 4, 8, 16, 32];
    minRange = [10^-3 10^-12];
    size_mcl = nan(length(Irange),1);
    for ii = 1:length(Irange)
        for jj = 1:length(Erange)
            for kk = 1:length(minRange)
                [~, clusts] = mymcl(intMatrix, Irange(ii), Erange(jj), minRange(kk), 20);
                nn = nan(size(clusts));
                for uu = 1:length(nn)
                    nn(uu) = length(clusts{uu});
                end
                nn = nn(nn>=3);
                size_mcl(ii,jj,kk) = median(nn);
            end
        end
    end

    
    % co+mcl
    pRange = [1, 50, 100, 500, 5000];
    densRange = [0, 0.1, 0.2, 0.3, 0.4];
    Irange = [1, 2, 4, 8, 16];
    minRange = [10^-3 10^-12];
    size_comcl = nan(length(Irange),1);
    for pp = 1:length(pRange)
        for ii = 1:length(densRange)
            for jj = 1:length(Irange)
                pars.best_p = pRange(pp);
                pars.best_dens = densRange(ii);
                pars.best_prec = 0.5;
                pars.best_I = Irange(jj);
                clusts = clustone_mcl(intMatrix, pars, java_path);
                nn = nan(size(clusts));
                for uu = 1:length(nn)
                    nn(uu) = length(clusts{uu});
                end
                nn = nn(nn>=3);
                size_comcl(pp,ii,jj) = median(nn);
            end
        end
    end

    
    save(fnsave, 'data', '-v7.3');
end



%% Write interactomes for functional analyses

% edit: don't bother with score, just write protA-protB

fnout = 'E:\Greg\ClusterReliable\data/interactomes_moredata_wscore.txt';
fid = fopen(fnout,'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\n','protA','protB','dataset','noise_mag','score');
cc = 0;
for ii = 1:size(data.chromnoise.interactome,1)
    this_dataset = ['us_' num2str(ii)];
    for jj = 1:size(data.chromnoise.interactome,2)
        this_noise_mag = data.chromnoiseRange(jj);
        for kk = 1:size(data.chromnoise.interactome{ii,jj})
            protA = data.Protein{ii}.NoIsoform{data.chromnoise.interactome{ii,jj}(kk,1)};
            protB = data.Protein{ii}.NoIsoform{data.chromnoise.interactome{ii,jj}(kk,2)};
            score = data.chromnoise.score{ii,jj}(data.chromnoise.interactome{ii,jj}(kk,1), data.chromnoise.interactome{ii,jj}(kk,2));
            fprintf(fid,'%s\t%s\t%s\t%6.4f\t%6.4f\n',...
                protA, protB, this_dataset, this_noise_mag, score);
        end
    end
end
fclose all;


% write shuffled corum
fnout = 'E:\Greg\ClusterReliable\data/interactomes_corum_shuffle.txt';
fid = fopen(fnout,'w');
fprintf(fid,'%s\t%s\t%s\t%s\n','protA','protB','dataset','noise_mag');
cc = 0;
for ii = 1:length(data.shufflenoiseRange)
    % ADD edges to corum-network
    intMatrix = shufflenetwork(data.corum.network, data.shufflenoiseRange(ii));
    [ia,ib] = find(intMatrix>0);
    I = ia>ib;
    ia = ia(I);
    ib = ib(I);
    for jj = 1:length(ia)
        protA = data.corum.proteins{ia(jj)};
        protB = data.corum.proteins{ib(jj)};
        fprintf(fid,'%s\t%s\t%s\t%6.4f\n',...
            protA, protB, 'corum', data.shufflenoiseRange(ii));
    end
end
fclose all;


% %interactome stored in data.chromnoise.score(:,1), i.e. no noise
% int2write.prots = cell(10^7,2);
% int2write.score = nan(10^7,1);
% int2write.dataset = nan(10^7,1);
% cc = 0;
% for ii = 1:size(data.chromnoise.score,1)
%     [ia,ib] = find(triu(data.chromnoise.score{ii,1}) > 0.5);
%     for jj = 1:length(ia)
%         cc = cc+1;
%         prots = sort({data.Protein{ii}.NoIsoform{ia(jj)} data.Protein{ii}.NoIsoform{ib(jj)}});
%         int2write.prots(cc,:) = prots;
%         int2write.score(cc) = data.chromnoise.score{ii,1}(ia(jj),ib(jj));
%         int2write.dataset(cc) = ii;
%     end
% end
% fnames = fieldnames(int2write);
% for ii = 1:length(fnames)
%     int2write.(fnames{ii}) = int2write.(fnames{ii})(1:cc,:);
% end
%
% fnout = 'E:\Greg\ClusterReliable\data/interactomes_moredata.txt';
% fid = fopen(fnout,'w');
% fprintf(fid,'%s\t%s\t%s\t%s\n','protA','protB','score','dataset');
% for ii = 1:cc
%     fprintf(fid,'%s\t%s\t%6.4f\t%d\n',...
%         int2write.prots{ii,1},int2write.prots{ii,2},...
%         int2write.score(ii),int2write.dataset(ii));
% end
% fclose all;


%% Write clusters for functional analysis (R)
% one cluster per row
% metadata:
%   - data type ('us' or 'corum')
%   - noise type ('chrom', 'network_add', 'network_remove', or 'network_shuffle')
%   - noise magnitude
%   - clustering algorithm

fnout = 'E:\Greg\ClusterReliable\data/clusters_wshuffle_moredata.txt';
fid = fopen(fnout, 'w');

fnames0 = {'add' 'remove' 'shuffle'};
fnames1 = {'co_mcl' 'co' 'mcl'};
fnames2 = {'mr' 'ga' 'sn' 'ppv' 'nmi' 'ari' 'coint' 'cocom'};

fprintf(fid,'%s\t%s\t%s\t%s\t%s\t\n','data_type','noise_type','noise_mag','algorithm','cluster');
% chromnoise
for ii = 1:length(fnames1) % ii = algorithm
    for jj = 1:size(data.chromnoise.(fnames1{ii}).cluster,1) % jj = dataset
        for kk = 1:size(data.chromnoise.(fnames1{ii}).cluster,2) % kk = noise mag
            these_clusters = data.chromnoise.(fnames1{ii}).cluster{jj,kk};
            for cc = 1:length(these_clusters)
                this_cluster = these_clusters{cc};
                if isempty(this_cluster); continue; end
                this_cluster = strjoin(data.Protein{jj}.MajorID_NoIsoforms(this_cluster), ';');
                fprintf(fid,'%s\t%s\t%6.4f\t%s\t%s\n',...
                    ['us_' num2str(jj)],'chrom',data.chromnoiseRange(kk),fnames1{ii},this_cluster);
            end
        end
    end
end
% us netnoise
for ii = 1:length(fnames1) % ii = algorithm
    for bb = 1:length(fnames0) % bb = add or remove
        for jj = 1:size(data.netnoise.(fnames0{bb}).(fnames1{ii}).cluster,1) % jj = dataset
            for kk = 1:size(data.netnoise.(fnames0{bb}).(fnames1{ii}).cluster,2) % kk = noise mag
                if ismember('add', fnames0{bb})
                    noisemag = data.addnoiseRange(kk);
                elseif ismember('remove', fnames0{bb})
                    noisemag = data.remnoiseRange(kk);
                elseif ismember('shuffle', fnames0{bb})
                    noisemag = data.shufflenoiseRange(kk);
                else
                    error('whaat')
                end
                these_clusters = data.netnoise.(fnames0{bb}).(fnames1{ii}).cluster{jj,kk};
                for cc = 1:length(these_clusters)
                    this_cluster = these_clusters{cc};
                    if isempty(this_cluster); continue; end
                    this_cluster = strjoin(data.Protein{jj}.MajorID_NoIsoforms(this_cluster), ';');
                    fprintf(fid,'%s\t%s\t%6.4f\t%s\t%s\n',...
                        ['us_' num2str(jj)],['network_' fnames0{bb}],noisemag,fnames1{ii},this_cluster);
                end
            end
        end
    end
end

% corum netnoise
for ii = 1:length(fnames1) % ii = algorithm
    for bb = 1:length(fnames0) % bb = add or remove
        for kk = 1:size(data.corum.(fnames0{bb}).(fnames1{ii}).cluster,2) % kk = noise mag
            if ismember('add', fnames0{bb})
                noisemag = data.addnoiseRange(kk);
            elseif ismember('remove', fnames0{bb})
                noisemag = data.remnoiseRange(kk);
            elseif ismember('shuffle', fnames0{bb})
                noisemag = data.shufflenoiseRange(kk);
            else
                error('whaat')
            end
            these_clusters = data.corum.(fnames0{bb}).(fnames1{ii}).cluster{kk};
            for cc = 1:length(these_clusters)
                this_cluster = these_clusters{cc};
                if isempty(this_cluster); continue; end
                this_cluster = strjoin(data.corum.proteins(this_cluster), ';');
                fprintf(fid,'%s\t%s\t%6.4f\t%s\t%s\n',...
                    'corum',['network_' fnames0{bb}],noisemag,fnames1{ii},this_cluster);
            end
        end
    end
end

fclose all;



% write shuffled_bad (not optimized parameters)

fnout = 'E:\Greg\ClusterReliable\data/clusters_wshuffle_notoptimal.txt';
fid = fopen(fnout, 'w');

fnames0 = {'shuffle_bad'};
fnames1 = {'co_mcl' 'co' 'mcl'};

fprintf(fid,'%s\t%s\t%s\t%s\t%s\n','data_type','noise_type','noise_mag','algorithm','cluster');
% us netnoise
for ii = 1:length(fnames1) % ii = algorithm
    for bb = 1:length(fnames0) % bb = add or remove
        for jj = 1:size(data.corum.(fnames0{bb}).(fnames1{ii}).cluster,1) % jj = dataset
            for kk = 1:size(data.corum.(fnames0{bb}).(fnames1{ii}).cluster,2) % kk = noise mag
                
                noisemag = data.shufflenoiseRange(kk);
                
                these_clusters = data.corum.(fnames0{bb}).(fnames1{ii}).cluster{jj,kk};
                for cc = 1:length(these_clusters)
                    this_cluster = these_clusters{cc};
                    if isempty(this_cluster); continue; end
                    this_cluster = strjoin(data.corum.proteins(this_cluster), ';');
                    fprintf(fid,'%s\t%s\t%6.4f\t%s\t%s\n',...
                        ['us_' num2str(jj)],['network_' fnames0{bb}],noisemag,fnames1{ii},this_cluster);
                end
            end
        end
    end
end



%% At the same noise level, how consistent are complexes?

% use CORUM, shuffle
clear data2
data2.noiseRange = [0 0.1 0.25 0.5];
algs = {'co_mcl' 'co' 'mcl'};
iterMax = 5;

for mm = 1:length(algs)
    data2.(algs{mm}).mr = cell(length(data2.noiseRange),1);%nan(iterMax, iterMax);
    data2.(algs{mm}).nmi = cell(length(data2.noiseRange),1);%nan(iterMax, iterMax);
    data2.(algs{mm}).ga = cell(length(data2.noiseRange),1);%nan(iterMax, iterMax);
    data2.(algs{mm}).sn = cell(length(data2.noiseRange),1);%nan(iterMax, iterMax);
    data2.(algs{mm}).ppv = cell(length(data2.noiseRange),1);%nan(iterMax, iterMax);
end

for ii = 1:length(data2.noiseRange)
    for jj = 1:iterMax
        disp(['shuffle ' num2str(data2.noiseRange(ii)) ' iter ' num2str(jj)])
        intMatrix = shufflenetwork(data.corum.network, data2.noiseRange(ii));
        
        % cluster with clusterone+mcl
        data2.co_mcl.cluster{ii,jj} = clustone_mcl(intMatrix, pars, java_path);
        
        % cluster with clusterone
        data2.co.cluster{ii,jj} = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path);
        
        % cluster with mcl
        [~, data2.mcl.cluster{ii,jj}] = mymcl(intMatrix, 1.8);
    end
    
    for mm = 1:length(algs)
        disp(['comparing networks ' algs{mm}])
        data2.(algs{mm}).mr{ii} = nan(iterMax, iterMax);
        data2.(algs{mm}).nmi{ii} = nan(iterMax, iterMax);
        data2.(algs{mm}).ga{ii} = nan(iterMax, iterMax);
        data2.(algs{mm}).sn{ii} = nan(iterMax, iterMax);
        data2.(algs{mm}).ppv{ii} = nan(iterMax, iterMax);
        if ismember(algs{mm},{'mcl'})
            data2.(algs{mm}).mr{ii} = nan(iterMax, iterMax);
        end
        for jj = 1:iterMax
            for kk = 1:iterMax
                % compare networks
                tmp = comparenetworks(data2.(algs{mm}).cluster{ii,jj},data2.(algs{mm}).cluster{ii,kk},...
                    intMatrix,data.corum.network);
                data2.(algs{mm}).mr{ii}(jj,kk) = tmp.mr;
                data2.(algs{mm}).nmi{ii}(jj,kk) = tmp.nmi;
                data2.(algs{mm}).ga{ii}(jj,kk) = tmp.ga;
                data2.(algs{mm}).sn{ii}(jj,kk) = tmp.sn;
                data2.(algs{mm}).ppv{ii}(jj,kk) = tmp.ppv;
                if ismember(algs{mm},{'mcl'})
                    data2.mcl.ari{ii}(jj,kk) = ari(data2.mcl.cluster{ii,jj},data2.mcl.cluster{ii,kk});
                end
            end
        end
    end
end



%% Make 100x chrom-noise iterations for a single-complex stability metric

iterMax = 100;
chromNoise = 0.1;

pars.best_p = 500;
pars.best_dens = 0.1;
pars.best_prec = 0.6;
pars.best_I = 4;

clear data2
data2.interactome = cell(length(data.file.name),iterMax);
data2.mcl.cluster = cell(length(data.file.name),iterMax);
data2.co.cluster = cell(length(data.file.name),iterMax);
data2.co_mcl.cluster = cell(length(data.file.name),iterMax);
fnames1 = {'co_mcl' 'co' 'mcl'};
tic
for ii = 1:length(data.file.name)
    for jj = 1:iterMax
        load(data.file.name{ii})
        disp(['Making score matrix... ' num2str(ii) ' ' num2str(jj)])
        % add noise to chroms
        if jj==1
            noised_chroms = data.Chromatograms{ii};
        else
            rnd = randn(size(data.Chromatograms{ii}));
            noised_chroms = data.Chromatograms{ii}.*exp(rnd*chromNoise);
        end
        
        % pre-process chromatograms
        chromclean = noised_chroms;
        I = find(isnan(chromclean));
        chromclean(I)= 0.05 * rand(size(I));
        
        % calc new Dist matrices
        clear tmpdist
        tmpdist.Euc = squareform(pdist(chromclean,'euclidean'));              % 1. Euclidean distance
        tmpdist.R = -1*corr(chromclean');                                     % 2. Cleaned chromatogram R^2
        [R,p] = corrcoef(chromclean','rows','pairwise');
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
        scoreMatrix = nan(size(TP_Matrix));
        for kk = 1:length(ia)
            scoreMatrix(ia(kk), ib(kk)) = prec(kk);
            scoreMatrix(ib(kk), ia(kk)) = prec(kk);
        end
        
        % get 50% interactome
        [ia,ib] = find(scoreMatrix >= 0.50);
        ibad = ia<ib;
        ia(ibad) = [];
        ib(ibad) = [];
        data2.interactome{ii,jj} = [ia ib];
        
        
        
        disp(['    chrom noise ' num2str(chromNoise*100) '%'])
        % make intMatrix
        intMatrix = scoreMatrix > 0.5;
        intMatrix(intMatrix<0) = 0;
        intMatrix = (intMatrix - nanmin(intMatrix(:))) / (max(intMatrix(:)) - min(intMatrix(:)));
        if nansum(intMatrix(:))==0; continue; end
        
        for kk = 1:length(fnames1)
            % cluster
            if ismember('co_mcl',fnames1(kk))
                disp('        co_mcl')
                data2.co_mcl.cluster{ii,jj} = clustone_mcl(intMatrix, pars, java_path);
            elseif ismember('co',fnames1(kk))
                disp('        co')
                data2.co.cluster{ii,jj} = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path);
            elseif ismember('mcl',fnames1(kk))
                disp('        mcl')
                [~,data2.mcl.cluster{ii,jj}] = mymcl(intMatrix, 1.8);
            end
        end
        toc
    end
end

% write data2 for analysis in R
fn = '../data/clusters_calcEA.txt';
fid = fopen(fn, 'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n','data_type','noise_type','iter','noise_mag','algorithm','cluster');
fns = {'mcl' 'co' 'co_mcl'};
for ii = 1:3 % clustering algorithm
    for jj = 1:8 % dataset
        for kk = 1:100 % iteration
            if kk==1
                noise_mag = 0;
            else
                noise_mag = 0.1;
            end
            clusters = data2.(fns{ii}).cluster{jj,kk};
            for mm = 1:length(clusters)
                this_cluster = data.Protein{jj}.NoIsoform(clusters{mm});
                this_cluster = [sprintf('%s;',this_cluster{1:end-1}),this_cluster{end}];
                fprintf(fid,'%d\t%s\t%d\t%6.4f\t%s\t%s\n',...
                    jj,'chrom',kk,noise_mag,fns{ii}, this_cluster);
            end
        end
    end
end
fclose all;

fn = '../data/interactomes_calcEA.txt';
fid = fopen(fn, 'w');
fprintf(fid,'%s\t%s\t%s\t%s\n','protA','protB','dataset','iter');
for jj = 1:8
    jj;
    for kk = 1:100
        prots = data.Protein{jj}.NoIsoform(data2.interactome{jj,kk});
        for ii = 1:length(prots)
            these_prots = sort(prots(ii,:));
            fprintf(fid,'%s\t%s\t%d\t%d\n',these_prots{1}, these_prots{2}, jj, kk);
        end
    end
end
fclose all;



%% Interactome FDR vs Ai analysis
% Make 25x corum clusters at a few interactome FDRs, e.g. FDR = [0,.1,.25]

iterMax = 10;

pars.best_p = 500;
pars.best_dens = 0.1;
pars.best_prec = 0.6;
pars.best_I = 4;

clear data3
data3.fdrRange = [0 0.1 0.25 0.5 .01 .02 .4 .75 .05 .15 .25];
data3.network = data.corum.network;
data3.proteins = data.corum.proteins;

fnames1 = {'co_mcl' 'co' 'mcl'};
for ii = 9:length(data3.fdrRange)
    disp(['fdr = ' num2str(data3.fdrRange(ii))])
    for jj = 1:iterMax
        % SHUFFLE edges of corum-network
        disp(['    iter = ' num2str(jj)])
        intMatrix = shufflenetwork(data.corum.network, data3.fdrRange(ii));
        
        disp('        clusterone+mcl')
        % cluster with clusterone+mcl
        data3.co_mcl.cluster{ii,jj} = clustone_mcl(intMatrix, pars, java_path);
        
        disp('        clusterone')
        % cluster with clusterone
        data3.co.cluster{ii,jj} = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path);
        
        disp('        mcl')
        % cluster with mcl
        [~, data3.mcl.cluster{ii,jj}] = mymcl(intMatrix, 1.8);
    end
end

% write clusters for analysis in R
% write data2 for analysis in R
fn = '../data/clusters_Ai_vs_fdr.txt';
fid = fopen(fn, 'w');
fprintf(fid,'%s\t%s\t%s\t%s\n','iter','noise_mag','algorithm','cluster');
fns = {'mcl' 'co' 'co_mcl'};
for kk = 1:size(data3.mcl.cluster,2) % iteration
    for jj = 1:length(data3.fdrRange)-1 % interactome fdr
        for ii = 1:3 % clustering algorithm
            clusters = data3.(fns{ii}).cluster{jj,kk};
            for mm = 1:length(clusters)
                this_cluster = data3.proteins(clusters{mm});
                this_cluster = [sprintf('%s;',this_cluster{1:end-1}),this_cluster{end}];
                fprintf(fid,'%d\t%6.4f\t%s\t%s\n',...
                    kk, data3.fdrRange(jj), fns{ii}, this_cluster);
            end
        end
    end
end
fclose all;