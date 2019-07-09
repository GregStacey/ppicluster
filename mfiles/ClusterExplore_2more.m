
% network noise / chromatogram noise
% dense/sparse

% TO DO:
%   1. In Craig netnoise, replace shufflenetwork with add+subtract.
%   1. In Craig, include all clustering algorithms.


%% 0. Get data

clear fns
fns{1} = 'E:\Greg\ClusterReliable\data/ChCh-Miner_durgbank-chem-chem.tsv'; % drug-drug network
fns{2} = 'E:\Greg\ClusterReliable\data/email-Eu-core.txt'; % email network
fn_gt = 'E:\Greg\ClusterReliable\data/email-Eu-core-department-labels.txt';

java_path = 'E:/Greg/ClusterReliable/java/cluster_one-1.0.jar';

fnsave = 'E:\Greg\ClusterReliable\data/data_facebook_save_01.mat';

try load(fnsave)
catch
    
    
    %% 1. turn networks into adjacency matrices
    
    clear data
    data.edges = cell(1,2);
    data.nodes = cell(1,2);
    data.adjacency = cell(1,2);
    seps = {'\t' ' '};
    for uu = 1:2
        
        % Read raw prot-chem file
        fid = fopen(fns{uu},'r');
        fgetl(fid); % header
        
        % read body of file
        cc = 0;
        data.edges{uu} = cell(10^6,2);
        while ~feof(fid)
            t1 = strsplit(fgetl(fid),seps{uu});
            cc = cc+1;
            data.edges{uu}{cc,1} = t1{1};
            data.edges{uu}{cc,2} = t1{2};
        end
        data.edges{uu} = data.edges{uu}(1:cc,:);
        fclose(fid);
        
        % All unique proteins
        data.nodes{uu} = unique(data.edges{uu}(1:cc, :));
        
        % Split complexes into pairwise list
        data.adjacency{uu} = zeros(length(data.nodes{uu}));
        [x,ia] = ismember(data.edges{uu}(:,1), data.nodes{uu});
        [x,ib] = ismember(data.edges{uu}(:,2), data.nodes{uu});
        for ii = 1:length(ia)
            data.adjacency{uu}(ia(ii),ib(ii)) = 1;
            data.adjacency{uu}(ib(ii),ia(ii)) = 1;
        end
        
        % remove self-interactions
        data.adjacency{uu}(eye(length(data.nodes{uu}))==1) = 0;
    end
    
    % email ground truth
    fid = fopen(fn_gt);
    email_ground = cell(42,1);
    while not(feof(fid))
        t1 = strsplit(fgetl(fid), ' ');
        node = str2double(t1{1});
        group = str2double(t1{2})+1;
        email_ground{group} = sort([email_ground{group} node]);
    end
    
    
    %% 2. optimize parameters
    
    x = -1;
    pars_co.p = [5000 1];
    pars_co.dens = [0.2 0.4];
    
    pars_mcl.I = [32 8];
    pars_mcl.E = [4 2];
    pars_mcl.minval = [10^-12 .001];
    
    pars_comcl.p = [100 1];
    pars_comcl.dens = [0.3 0.1];
    pars_comcl.I = [1 8];
    
    if 0
        uu = 1;
        intMatrix = data.adjacency{uu};
        
        % co
        pRange = [1, 50, 100, 500, 5000];
        densRange = [0, 0.1, 0.2, 0.3, 0.4];
        mr_co = nan(length(pRange), length(densRange));
        for ii = 1:length(pRange)
            for jj = 1:length(densRange)
                clusts = clusterone_java(intMatrix, pRange(ii), densRange(jj), java_path);
                if uu==2 % use mr with a ground truth
                    mr_co(ii,jj) = matchingratio(clusts, email_ground);
                elseif uu==1 % use silhouette without a ground truth
                    mr_co(ii,jj) = mysilhouette(clusts, intMatrix);
                end
            end
        end
        I = find(mr_co == nanmax(mr_co(:)));
        [ia,ib] = ind2sub(size(mr_co), I);
        pars_co.p = pRange(ia);
        pars_co.dens = densRange(ib);
        
        % mcl
        Irange = [2, 4, 8, 16, 32];
        Erange = [2, 4, 8, 16, 32];
        minRange = [10^-3 10^-12];
        mr_mcl = nan(length(Irange),1);
        nn_mcl = nan(length(Irange),1);
        for ii = 1:length(Irange)
            for jj = 1:length(Erange)
                for kk = 1:length(minRange)
                    [~, clusts] = mymcl(intMatrix, Irange(ii), Erange(jj), minRange(kk), 20);
                    nn_mcl(ii,jj,kk) = length(clusts);
                    if uu==2
                        mr_mcl(ii,jj,kk) = matchingratio(clusts, email_ground);
                    elseif uu==1
                        mr_mcl(ii,jj,kk) = mysilhouette(clusts, intMatrix);
                    end
                end
            end
        end
        I = find(mr_mcl == nanmax(mr_mcl(:)));
        [ia,ib,ic] = ind2sub(size(mr_mcl), I);
        pars_mcl.I = Irange(ia);
        pars_mcl.E = Erange(ib);
        pars_mcl.minval = minRange(ic);
        
        % co+mcl
        pRange = [1, 50, 100, 500, 5000];
        densRange = [0, 0.1, 0.2, 0.3, 0.4];
        Irange = [1, 2, 4, 8, 16];
        minRange = [10^-3 10^-12];
        mr_comcl = nan(length(Irange),1);
        for pp = 1:length(pRange)
            for ii = 1:length(densRange)
                for jj = 1:length(Irange)
                    pars.best_p = pRange(pp);
                    pars.best_dens = densRange(ii);
                    pars.best_prec = 0.5;
                    pars.best_I = Irange(jj);
                    clusts = clustone_mcl(intMatrix, pars, java_path);
                    if uu==2
                        mr_comcl(pp,ii,jj) = matchingratio(clusts, email_ground);
                    elseif uu==1
                        mr_comcl(pp,ii,jj) = mysilhouette(clusts, intMatrix);
                    end
                end
            end
        end
        I = find(mr_comcl == nanmax(mr_comcl(:)));
        [ia,ib,ic] = ind2sub(size(mr_comcl), I);
        pars_comcl.p = pRange(ia);
        pars_comcl.dens = densRange(ib);
        pars_comcl.I = Irange(ic);
    end
    
    
    %% 2. Cluster networks

    data.shufflenoiseRange = [0 0.01 0.02 0.05 .1 0.15 .25 .5 1];
    data.clusters_mcl = cell(length(data.shufflenoiseRange), 2);
    data.clusters_co = cell(length(data.shufflenoiseRange), 2);
    data.clusters_comcl = cell(length(data.shufflenoiseRange), 2);
    
    for uu = 1:2
        for ii = 1:length(data.shufflenoiseRange)
            disp('Clustering corum... ')
            
            % SHUFFLE edges of corum-network
            disp(['    shuffle ' num2str(data.shufflenoiseRange(ii))])
            intMatrix = shufflenetwork(data.adjacency{uu}, data.shufflenoiseRange(ii));
            
            disp('        clusterone+mcl')
            % cluster with clusterone+mcl
            pars.best_p = pars_comcl.p(uu);
            pars.best_dens =  pars_comcl.dens(uu);
            pars.best_prec = 0.1;
            pars.best_I =  pars_comcl.I(uu);
            data.clusters_comcl{ii, uu} = clustone_mcl(intMatrix, pars, java_path);
            
            disp('        clusterone')
            % cluster with clusterone
            data.clusters_co{ii, uu} = clusterone_java(intMatrix, pars_co.p(uu), pars_co.dens(uu), java_path);
            
            disp('        mcl')
            % cluster with mcl
            %[~, data.clusters_mcl{ii, uu}] = mymcl(intMatrix, pars_mcl.I(uu), pars_mcl.E(uu), pars_mcl.minval(uu), 20);
            if ii>5
            [~, data.clusters_mcl{ii, uu}] = mymcl(intMatrix, 8, 4, pars_mcl.minval(uu), 20);
            else
            [~, data.clusters_mcl{ii, uu}] = mymcl(intMatrix, 8, 2, pars_mcl.minval(uu), 20);
            end
        end
    end
    
    save(fnsave, 'data', '-v7.3');
end


%% Write clusters

fnames = {'clusters_mcl' 'clusters_comcl' 'clusters_co'};
netnames = {'chem' 'email'};
algnames = {'mcl' 'co_mcl' 'co'};

fnout = 'E:\Greg\ClusterReliable\data/clusters_chem_email.txt';
fid = fopen(fnout, 'w');
fprintf(fid,'%s\t%s\t%s\t%s\n','network','noise_mag','algorithm','cluster');

for ii = 1:length(fnames) % ii = algorithm
    for uu = 1:2 % network
        for kk = 1:length(data.shufflenoiseRange) % kk = noise mag
            
            noisemag = data.shufflenoiseRange(kk);
            
            these_clusters = data.(fnames{ii}){kk,uu};
            for cc = 1:length(these_clusters)
                this_cluster = these_clusters{cc};
                if isempty(this_cluster); continue; end
                this_cluster = strjoin(data.nodes{uu}(this_cluster), ';');
                fprintf(fid,'%s\t%6.4f\t%s\t%s\n',...
                    netnames{uu},noisemag,algnames{ii},this_cluster);
            end
        end
    end
end
fclose all;

