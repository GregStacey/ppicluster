
% network noise / chromatogram noise
% dense/sparse

% clear data
% data.file.name = data.file;
% data.file.rep = data.rep;
% data.file.channel = data.channel;
% data.Chromatograms = data.Chromatograms;
% data.Protein = data.Protein;
% data.chromnoise.Chromatograms = data.Chrom_noise;
% data.chromnoise.interactome = data.chromnoise.interactome;
% data.chromnoise.sparse.cluster = data.cluster_chrom_noise;
% data.chromnoise.sparse.mr = data.mr_chrom_noise;
% data.chromnoise.sparse.ga = data.ga_chrom_noise;
% data.chromnoise.sparse.coint = data.cointprob_chrom_noise;
% data.chromnoise.sparse.cocom = data.cocomprob_chrom_noise;
% data.chromnoise.dense.cluster = data.cluster_chrom_noise_dense;
% data.chromnoise.dense.mr = data.mr_chrom_noise_dense;
% data.chromnoise.dense.ga = data.ga_chrom_noise_dense;
% data.chromnoise.dense.coint = data.cointprob_chrom_noise_dense;
% data.chromnoise.dense.cocom = data.cocomprob_chrom_noise_dense;
% data.netnoise.dense.cluster = data.cluster_network_noise;
% data.netnoise.dense.mr = data.mr_network_noise;
% data.netnoise.dense.ga = data.ga_network_noise;
% data.netnoise.dense.coint = data.cointprob_network_noise;
% data.netnoise.dense.cocom = data.cocomprob_network_noise;
% data.netnoise.sparse.cluster = data.cluster_network_noise_dense;
% data.netnoise.sparse.mr = data.mr_network_noise_dense;
% data.netnoise.sparse.ga = data.ga_network_noise_dense;
% data.netnoise.sparse.coint = data.cointprob_network_noise_dense;
% data.netnoise.sparse.cocom = data.cocomprob_network_noise_dense;
% data = data;


%% 0. Get data

datadir = 'E:/Greg/PCP_SILAC/Runs/Craig/Final_interactome_13k/Output/tmp/';
datanames = {'Craig'};
tmp = dir([datadir 'data*rep*chan*.mat']);

fnsave = 'E:\Greg\ClusterReliable\data/data_save_05.mat';

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
    
    
    
    %% 2a. Cluster (chrom noise, sparse unweighted)
    
    pars.best_p = 500;
    pars.best_dens = 0.1;
    pars.best_prec = 0.6;
    pars.best_I = 4;
    
    data.chromnoise.sparse.cluster = cell(size(data.chromnoise.score));
    data.chromnoise.sparse.mr = nan(size(data.chromnoise.score));
    data.chromnoise.sparse.ga = nan(size(data.chromnoise.score));
    data.chromnoise.sparse.cocom = nan(size(data.chromnoise.score));
    data.chromnoise.sparse.coint = nan(size(data.chromnoise.score));
    for ii = 1:length(data.file.name)
        % chrom_noise clusters
        for jj = 1:size(data.chromnoise.score,2)
            disp(['Clustering chrom_noise... ' num2str(ii) ' ' num2str(jj)])
            % make intMatrix
            intMatrix = data.chromnoise.score{ii,jj} > 0.5;
            intMatrix(intMatrix<0) = 0;
            intMatrix = (intMatrix - nanmin(intMatrix(:))) / (nanmax(intMatrix(:)) - nanmin(intMatrix(:)));
            
            % cluster
            data.chromnoise.sparse.cluster{ii,jj} = clustone_mcl(intMatrix, pars);
            
            disp(['Calculating MR... ' num2str(ii) ' ' num2str(jj)])
            
            if isempty(data.chromnoise.sparse.cluster{ii,1}); continue; end
            if isempty(data.chromnoise.sparse.cluster{ii,jj}); continue; end
            
            % compare networks
            tmp = comparenetworks(data.chromnoise.sparse.cluster{ii,jj},...
                data.chromnoise.sparse.cluster{ii,1},...
                data.chromnoise.score{ii,jj},data.chromnoise.score{ii,1});
            data.chromnoise.sparse.mr(ii,jj) = tmp.mr;
            data.chromnoise.sparse.ga(ii,jj) = tmp.ga;
            data.chromnoise.sparse.coint(ii,jj) = tmp.coint;
            data.chromnoise.sparse.cocom(ii,jj) = tmp.cocom;
        end
    end
    
    
    
    %% 2b. Cluster (chrom_noise, dense weighted)
    pars.best_p = 50000;
    pars.best_dens = 0;
    pars.best_prec = 0;
    pars.best_I = 100;
    
    data.chromnoise.dense.cluster = cell(size(data.chromnoise.score));
    data.chromnoise.dense.mr = nan(size(data.chromnoise.score));
    data.chromnoise.dense.ga = nan(size(data.chromnoise.score));
    data.chromnoise.dense.cocom = nan(size(data.chromnoise.score));
    data.chromnoise.dense.coint = nan(size(data.chromnoise.score));
    for ii = 1:length(data.file.name)
        % chrom_noise clusters
        for jj = 1:size(data.chromnoise.score,2)
            disp(['Clustering chrom_noise... ' num2str(ii) ' ' num2str(jj) ' (dense)'])
            % make intMatrix
            intMatrix = data.chromnoise.score{ii,jj} .^ 2;
            intMatrix(isnan(intMatrix)) = 0;
            
            % cluster
            data.chromnoise.dense.cluster{ii,jj} = clustone_mcl(intMatrix, pars);
            
            disp(['Calculating MR... ' num2str(ii) ' ' num2str(jj) ' (dense)'])
            
            if isempty(data.chromnoise.dense.cluster{ii,1}); continue; end
            if isempty(data.chromnoise.dense.cluster{ii,jj}); continue; end
            
            % compare networks
            tmp = comparenetworks(data.chromnoise.dense.cluster{ii,jj},...
                data.chromnoise.dense.cluster{ii,1},...
                intMatrix,data.chromnoise.score{ii,1});
            data.chromnoise.dense.mr(ii,jj) = tmp.mr;
            data.chromnoise.dense.ga(ii,jj) = tmp.ga;
            data.chromnoise.dense.coint(ii,jj) = tmp.coint;
            data.chromnoise.dense.cocom(ii,jj) = tmp.cocom;
        end
    end
    
    
    %% 2c. Cluster (network noise, sparse unweighted)
    
    pars.best_p = 500;
    pars.best_dens = 0.1;
    pars.best_prec = 0.6;
    pars.best_I = 4;
    
    data.netnoise.sparse.cluster = cell(size(data.chromnoise.score));
    data.netnoise.sparse.mr = nan(size(data.chromnoise.score));
    data.netnoise.sparse.ga = nan(size(data.chromnoise.score));
    data.netnoise.sparse.cocom = nan(size(data.chromnoise.score));
    data.netnoise.sparse.coint = nan(size(data.chromnoise.score));
    for ii = 1:length(data.file.name)
        % chrom_noise clusters
        for jj = 1:size(data.chromnoise.score,2)
            disp(['Clustering chrom_noise... ' num2str(ii) ' ' num2str(jj)])
            % make intMatrix
            intMatrix = data.chromnoise.score{ii,1};
            intMatrix = addnoise2network(intMatrix, mRange(jj));
            intMatrix = intMatrix > 0.5;
            intMatrix(intMatrix<0) = 0;
            intMatrix = (intMatrix - nanmin(intMatrix(:))) / (nanmax(intMatrix(:)) - nanmin(intMatrix(:)));
            
            % cluster
            data.netnoise.sparse.cluster{ii,jj} = clustone_mcl(intMatrix, pars);
            
            disp(['Calculating MR... ' num2str(ii) ' ' num2str(jj)])
            
            if isempty(data.netnoise.sparse.cluster{ii,1}); continue; end
            if isempty(data.netnoise.sparse.cluster{ii,jj}); continue; end
            
            % compare networks
            tmp = comparenetworks(data.netnoise.sparse.cluster{ii,jj},...
                data.netnoise.sparse.cluster{ii,1},...
                intMatrix,data.chromnoise.score{ii,1});
            data.netnoise.sparse.mr(ii,jj) = tmp.mr;
            data.netnoise.sparse.ga(ii,jj) = tmp.ga;
            data.netnoise.sparse.coint(ii,jj) = tmp.coint;
            data.netnoise.sparse.cocom(ii,jj) = tmp.cocom;
        end
    end
    
    
    %% 2d. Cluster (network noise, dense)
    pars.best_p = 50000;
    pars.best_dens = 0;
    pars.best_prec = 0;
    pars.best_I = 100;
    
    data.netnoise.dense.cluster = cell(size(data.chromnoise.score));
    data.netnoise.dense.mr = nan(size(data.chromnoise.score));
    data.netnoise.dense.ga = nan(size(data.chromnoise.score));
    data.netnoise.dense.cocom = nan(size(data.chromnoise.score));
    data.netnoise.dense.coint = nan(size(data.chromnoise.score));
    for ii = 1:length(data.file.name)
        % network_noise clusters
        for jj = 1:size(data.chromnoise.score,2)
            disp(['Clustering network_noise... ' num2str(ii) ' ' num2str(jj) ' (dense)'])
            % make intMatrix
            intMatrix = data.chromnoise.score{ii,1} .^ 2;
            intMatrix = addnoise2network(intMatrix, mRange(jj));
            intMatrix(isnan(intMatrix)) = 0;
            
            % cluster
            data.netnoise.dense.cluster{ii,jj} = clustone_mcl(intMatrix, pars);
            
            disp(['Calculating MR... ' num2str(ii) ' ' num2str(jj) ' (dense)'])
            
            if isempty(data.netnoise.dense.cluster{ii,1}); continue; end
            if isempty(data.netnoise.dense.cluster{ii,jj}); continue; end
            
            % compare networks
            tmp = comparenetworks(data.netnoise.dense.cluster{ii,jj},...
                data.netnoise.dense.cluster{ii,1},...
                intMatrix,data.chromnoise.score{ii,1} .^ 2);
            data.netnoise.dense.mr(ii,jj) = tmp.mr;
            data.netnoise.dense.ga(ii,jj) = tmp.ga;
            data.netnoise.dense.coint(ii,jj) = tmp.coint;
            data.netnoise.dense.cocom(ii,jj) = tmp.cocom;
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
        data.corum.add.co_mcl.cluster{ii} = clustone_mcl(intMatrix, pars);
        % compare networks
        tmp = comparenetworks(data.corum.add.co_mcl.cluster{ii},data.corum.add.co_mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.add.co_mcl.mr(ii) = tmp.mr;
        data.corum.add.co_mcl.ga(ii) = tmp.ga;
        data.corum.add.co_mcl.coint(ii) = tmp.coint;
        data.corum.add.co_mcl.cocom(ii) = tmp.cocom;
        
        disp('        clusterone')
        % cluster with clusterone
        data.corum.add.co.cluster{ii} = myclusterone(intMatrix, pars.best_p, 0);
        % compare networks
        tmp = comparenetworks(data.corum.add.co.cluster{ii},data.corum.add.co.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.add.co.mr(ii) = tmp.mr;
        data.corum.add.co.ga(ii) = tmp.ga;
        data.corum.add.co.coint(ii) = tmp.coint;
        data.corum.add.co.cocom(ii) = tmp.cocom;
        
        disp('        mcl')
        % cluster with mcl
        [~, data.corum.add.mcl.cluster{ii}] = mymcl(intMatrix, 1.8);
        % compare networks
        tmp = comparenetworks(data.corum.add.mcl.cluster{ii},data.corum.add.mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.add.mcl.mr(ii) = tmp.mr;
        data.corum.add.mcl.ga(ii) = tmp.ga;
        data.corum.add.mcl.coint(ii) = tmp.coint;
        data.corum.add.mcl.cocom(ii) = tmp.cocom;
        
        
        % REMOVE edges to corum-network
        disp(['    remove ' num2str(cRange(ii))])
        intMatrix = addremovenetwork(data.corum.network, -1 * cRange(ii));
        
        disp('        clusterone+mcl')
        % cluster with clusterone+mcl
        data.corum.remove.co_mcl.cluster{ii} = clustone_mcl(intMatrix, pars);
        % compare networks
        tmp = comparenetworks(data.corum.remove.co_mcl.cluster{ii},data.corum.remove.co_mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.remove.co_mcl.mr(ii) = tmp.mr;
        data.corum.remove.co_mcl.ga(ii) = tmp.ga;
        data.corum.remove.co_mcl.coint(ii) = tmp.coint;
        data.corum.remove.co_mcl.cocom(ii) = tmp.cocom;
        
        disp('        clusterone')
        % cluster with clusterone
        data.corum.remove.co.cluster{ii} = myclusterone(intMatrix, pars.best_p, 0);
        % compare networks
        tmp = comparenetworks(data.corum.remove.co.cluster{ii},data.corum.remove.co.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.remove.co.mr(ii) = tmp.mr;
        data.corum.remove.co.ga(ii) = tmp.ga;
        data.corum.remove.co.coint(ii) = tmp.coint;
        data.corum.remove.co.cocom(ii) = tmp.cocom;
        
        disp('        mcl')
        % cluster with mcl
        [~, data.corum.remove.mcl.cluster{ii}] = mymcl(intMatrix, 1.8);
        % compare networks
        tmp = comparenetworks(data.corum.remove.mcl.cluster{ii},data.corum.remove.mcl.cluster{1},...
            intMatrix,data.corum.network);
        data.corum.remove.mcl.mr(ii) = tmp.mr;
        data.corum.remove.mcl.ga(ii) = tmp.ga;
        data.corum.remove.mcl.coint(ii) = tmp.coint;
        data.corum.remove.mcl.cocom(ii) = tmp.cocom;
    end
    
    save(fnsave, 'data');
end



