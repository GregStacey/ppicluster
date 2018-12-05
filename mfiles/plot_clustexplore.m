%% initialize

figdir = 'E:/Greg/ClusterReliable\figures/';

% colormap - jet with grey at the bottom
cmap = colormap('jet');
cmap = [.65 .65 .65; cmap];
colormap(cmap)


%% Plot

sRange = [0 0.01 0.02 0.05 .1 0.15 .25 .5 1];


figure,hold on
plot(data.chromnoise.Chromatograms{1,1}(395,:),'k')
plot(data.chromnoise.Chromatograms{1,5}(395,:),'b')
plot(data.chromnoise.Chromatograms{1,7}(395,:),'r')
xlabel('Fraction')
ylabel('SILAC ratio')
legend('Original','10% noise','25% noise')
fn = 'E:\Greg\ClusterReliable\figures/chromNoise.png';
saveas(gcf,fn)


figure,subplot(2,1,1)
plot(sRange,data.chromnoise.sparse.coint([2 3 4 6 7 8],:)')
xlim([0 .5])
ylabel('Co-interactome probability')
title('Chrom noise, co-int vs co-com')
subplot(2,1,2)
plot(sRange,data.chromnoise.sparse.cocom([2 3 4 6 7 8],:)')
xlim([0 .5])
ylabel('Co-complex probability')
xlabel('Noise added (% amplitude)')
fn = 'E:\Greg\ClusterReliable\figures/chromNoise_coint_cocomp.png';
saveas(gcf,fn)


% Chrom noise, Dense vs sparse
I = [2 3 4 6 7 8];
figure,subplot(1,3,1), hold on
plot(sRange,data.chromnoise.sparse.cocom(I,:)','k')
plot(sRange,data.chromnoise.dense.cocom(I,:)', 'r')
axis([0 .5 0 1])
ylabel('Co-complex probability')
xlabel('Noise added (% amplitude)')
title('Chrom noise, dense vs sparse')

subplot(1,3,2), hold on
plot(sRange,data.chromnoise.sparse.ga(I,:)', 'k')
plot(sRange,data.chromnoise.dense.ga(I,:)', 'r')
axis([0 .5 0 1])
ylabel('Geometric accuracy')
xlabel('Noise added (% amplitude)')

subplot(1,3,3), hold on
plot(0,0,'k')
plot(0,0,'r')
plot(sRange,data.chromnoise.sparse.mr(I,:)', 'k')
plot(sRange,data.chromnoise.dense.mr(I,:)', 'r')
axis([0 .5 0 1])
ylabel('Matching ratio')
xlabel('Noise added (% amplitude)')
legend('Sparse','Dense')
set(gcf,'units','normalized','paperunits','normalized',...
    'position',[.1 .1 .6 .3],'paperposition',[.1 .1 .6 .3])
fn = 'E:\Greg\ClusterReliable\figures/chromNoise_dense_vs_sparse.png';
saveas(gcf,fn)


% Co-int vs Co-comp (clustering amplifies noise)
figure, hold on
scatter(data.chromnoise.sparse.coint(:), data.chromnoise.sparse.cocom(:), 40, 'filled')
grid on
plot([0 1], [0 1], '--k')
xlabel('Co-interactome probability')
ylabel('Co-complex probability')


% Network noise, Dense vs sparse
%I = [2 3 4 6 7 8];
I = 1:8;
figure,subplot(1,3,1), hold on
plot(mRange,data.netnoise.dense.cocom(I,:)','k')
plot(mRange,data.netnoise.sparse.cocom(I,:)', 'r')
axis([0 .5 0 1])
ylabel('Co-complex probability')
xlabel('Noise added (% shuffled)')
title('Network noise, dense vs sparse')
subplot(1,3,2), hold on
plot(mRange,data.netnoise.dense.ga(I,:)', 'k')
plot(mRange,data.netnoise.sparse.ga(I,:)', 'r')
axis([0 .5 0 1])
ylabel('Geometric accuracy')
xlabel('Noise added (% shuffled)')
subplot(1,3,3), hold on
plot(0,0,'k')
plot(0,0,'r')
plot(mRange,data.netnoise.dense.mr(I,:)', 'k')
plot(mRange,data.netnoise.sparse.mr(I,:)', 'r')
axis([0 .5 0 1])
ylabel('Matching ratio')
xlabel('Noise added (% shuffled)')
legend('Sparse','Dense')
set(gcf,'units','normalized','paperunits','normalized',...
    'position',[.1 .1 .6 .3],'paperposition',[.1 .1 .6 .3])



% Corum-noise
figure, subplot(2,4,1), hold on
plot(cRange, data.corum.add.co_mcl.mr,'k')
plot(cRange, data.corum.add.co.mr,'r')
plot(cRange, data.corum.add.mcl.mr,'b')
legend('co+mcl','co','mcl')
title('MMR')

subplot(2,4,2), hold on
plot(cRange, data.corum.add.co_mcl.ga,'k')
plot(cRange, data.corum.add.co.ga,'r')
plot(cRange, data.corum.add.mcl.ga,'b')
title('GA')

subplot(2,4,3), hold on
plot(cRange, data.corum.add.co_mcl.coint,'k')
plot(cRange, data.corum.add.co.coint,'r')
plot(cRange, data.corum.add.mcl.coint,'b')
title('Co-interactome prob')

subplot(2,4,4), hold on
plot(cRange, data.corum.add.co_mcl.cocom,'k')
plot(cRange, data.corum.add.co.cocom,'r')
plot(cRange, data.corum.add.mcl.cocom,'b')
title('Co-complex prob')

subplot(2,4,5), hold on
plot(cRange(1:8), data.corum.remove.co_mcl.mr,'k')
plot(cRange(1:8), data.corum.remove.co.mr,'r')
plot(cRange, data.corum.remove.mcl.mr,'b')

subplot(2,4,6), hold on
plot(cRange(1:8), data.corum.remove.co_mcl.ga,'k')
plot(cRange(1:8), data.corum.remove.co.ga,'r')
plot(cRange, data.corum.remove.mcl.ga,'b')

subplot(2,4,7), hold on
plot(cRange(1:8), data.corum.remove.co_mcl.coint,'k')
plot(cRange(1:8), data.corum.remove.co.coint,'r')
plot(cRange, data.corum.remove.mcl.coint,'b')

subplot(2,4,8), hold on
plot(cRange(1:8), data.corum.remove.co_mcl.cocom,'k')
plot(cRange(1:8), data.corum.remove.co.cocom,'r')
plot(cRange, data.corum.remove.mcl.cocom,'b')




%% Illustrate clustering with an "animation"

% example complexes: ribosome, proteasome, CCT

fn = 'E:/Greg/ClusterReliable/data/allComplexes.txt';
fid = fopen(fn);
fgetl(fid);
corum = cell(10000, 3);
Nmembers = nan(10000,1);
cc = 0;
while not(feof(fid))
    cc = cc+1;
    t1 = strsplit(fgetl(fid), '\t');
    corum{cc,1} = t1{2}; % complex name
    corum{cc,2} = strsplit(t1{6},';'); % members, split
    corum{cc,3} = t1{3}; % organism
    Nmembers(cc) = length(corum{cc,2});
end
corum = corum(1:cc,:);
Nmembers = Nmembers(1:cc);
Nmembers = Nmembers(ismember(corum(:,3), 'Human'));
corum = corum(ismember(corum(:,3),'Human'),:);
[x,unqprots] = corum2network(fn);


I0 = randsample(find(Nmembers>5), 10);
prots = [];
for ii = 1:length(I0)
    prots = [prots corum{I0(ii),2}];
end

Iprots = nan(size(prots));
for ii = 1:length(Iprots)
    Iprots(ii) = find(contains(unqprots, prots{ii}));
end


% ground truth
M = data.corum.network(Iprots, Iprots);
M_cluster = zeros(size(M))-1;
saved_color_complex = nan(length(I0),1);
for ii = 1:length(Iprots)
    prota = unqprots{Iprots(ii)};
    I = cellfun(@(x) ismember(x, prota), corum(I0,2), 'UniformOutput', 0);
    I = cellfun(@(x) sum(x)>0, I);
    this_complex_a = find(I);
    for jj = 1:length(Iprots)
        if ii>=jj; continue; end
        protb = unqprots{Iprots(jj)};
        I = cellfun(@(x) ismember(x, protb), corum(I0,2), 'UniformOutput', 0);
        I = cellfun(@(x) sum(x)>0, I);
        this_complex_b = find(I);
        color_complex = intersect(this_complex_a, this_complex_b);
        if not(isempty(color_complex))
            saved_color_complex(jj) = color_complex(end);
            M_cluster(ii,jj) = color_complex(end);
            M_cluster(jj,ii) = color_complex(end);
        end
    end
end
% ground truth network
figure
imagesc(M)
axis xy square
colormap bone
set(gcf,'paperunits','centimeters','paperposition',[1 1 5 6],...
    'units','centimeters','position',[1 1 5 6])
sf = [figdir 'figure01_netw0_v01.png'];
saveas(gcf,sf);
% ground truth complexes
figure
imagesc(M+M_cluster)
colormap(cmap)
axis xy square
set(gcf,'paperunits','centimeters','paperposition',[1 1 5 6],...
    'units','centimeters','position',[1 1 5 6])
sf = [figdir 'figure01_clust0_v01.png'];
saveas(gcf,sf);


for ii = 1:8
    C = data.corum.add.mcl.cluster{ii};
    I = cellfun(@(x) ismember(x, Iprots), C, 'UniformOutput', 0);
    I = cellfun(@(x) sum(x)>0, I);
    goodClusters = C(I);
    
    % match the colours of goodClusters to corum{I0,2}
    JJ = nan(length(goodClusters), length(I0));
    for jj = 1:length(goodClusters)
        this_jj = unqprots(goodClusters{jj});
        for kk = 1:length(I0)
            this_kk = corum{I0(kk),2};
            JJ(kk,jj) = length(intersect(this_jj, this_kk));% / ...
                %length(unique([this_kk this_jj']));
        end
    end
    % resolve multiple predicted mathcing to the same reference
    Ibad = find(nansum(JJ>0,2)>1);
    for jj = 1:length(Ibad)
        mx = nanmax(JJ(Ibad(jj),:));
        I = JJ(Ibad(jj),:) == mx;
        JJ(Ibad(jj), not(I)) = nan;
    end
    [x,goodClusters_color] = max(JJ);
    % resolve complexes not matched to a reference
    Ibad = find(nansum(JJ)==0);
    for jj = 1:length(Ibad)
        goodClusters_color(Ibad(jj)) = min(goodClusters_color) + (jj-1) + 0.5;
    end
    
    cluster_color = zeros(length(Iprots), length(Iprots))-2;
    for jj = 1:length(goodClusters)
        this_cluster = intersect(goodClusters{jj}, Iprots);
        for kk = 1:length(this_cluster)
            for mm = 1:length(this_cluster)
                if kk>=mm; continue; end
                ia = find(Iprots==this_cluster(kk));
                ib = find(Iprots==this_cluster(mm));
                cluster_color(ia,ib) = goodClusters_color(jj);
                cluster_color(ib,ia) = goodClusters_color(jj);
            end
        end
    end
    figure
    M = addremovenetwork(data.corum.network(Iprots,Iprots), 1 * cRange(ii));
    imagesc(M)
    colormap('bone')
    axis xy square
    set(gcf,'paperunits','centimeters','paperposition',[1 1 5 6],...
        'units','centimeters','position',[1 1 5 6])
    sf = [figdir 'figure01_netw' num2str(ii) '_v01.png'];
    saveas(gcf,sf);
    
    figure
    imagesc(cluster_color+1)
    colormap(cmap)
    axis xy square
    set(gcf,'paperunits','centimeters','paperposition',[1 1 5 6],...
        'units','centimeters','position',[1 1 5 6])
    sf = [figdir 'figure01_mcl' num2str(ii) '_v01.png'];
    saveas(gcf,sf);
end


%% Write intMatrix for confirmation by clusterone.java

clear MM
% no noise
MM{1} = data.corum.network(Iprots, Iprots);
% 10% noise
MM{2} = addremovenetwork(data.corum.network(Iprots,Iprots), 0.1);
% 100% noise
MM{3} = addremovenetwork(data.corum.network(Iprots,Iprots), 0.2);
fns = {'0' '10' '100'};

for ii = 1:3
    fnout = ['E:\Greg\ClusterReliable\data/network4cojavavalidation_' num2str(ii) '.txt'];
    [ia,ib] = find(MM{ii}>0);
    fid = fopen(fnout,'w');
    for jj = 1:length(ia)
        fprintf(fid,'%d %d 1\n',ia(jj),ib(jj));
    end
    fclose(fid);
end

