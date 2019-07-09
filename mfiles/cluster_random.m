% cheeky!
% cluster random data

for iter = 1:25
    NN = ceil(100 + rand*500);
    pars.best_p = 250 + round(rand*500);
    pars.best_dens = 0.1;
    pars.best_prec = 0.2;
    pars.best_I = 2 + round(rand*5);
    java_path = 'E:/Greg/ClusterReliable/java/cluster_one-1.0.jar';
    
    
    % get N random corum proteins
    fn = '../data/allComplexes.txt';
    fid = fopen(fn);
    fgetl(fid);
    corum = cell(3000,1);
    cc = 0;
    while not(feof(fid))
        t1 = strsplit(fgetl(fid), '\t');
        cc = cc+1;
        corum{cc} = strsplit(t1{6}, ';');
    end
    corum = corum(1:cc);
    fclose all;
    
    prots = unique([corum{:}]);
    these_prots = randsample(prots, NN);
    
    
    % give each a random chromatogram
    chroms = rand(length(these_prots), 50);
    
    % cluster random data
    intMatrix = corr(chroms').^2;
    clusters = clusterone_java(intMatrix, pars.best_p, pars.best_dens, java_path);
    Iclusters = zeros(NN,1);
    for ii = 1:length(clusters)
        Iclusters(clusters{ii}) = ii;
    end
    length(unique([clusters{:}]))
    
    % export for GO enrichment in R
    fn = ['../data/random_clusters' num2str(iter) '.txt'];
    fid = fopen(fn, 'w');
    fprintf(fid, 'protein\tcluster');
    for ii = 1:size(chroms,2)
        fprintf(fid,'\t%s', ['f' num2str(ii)]);
    end
    fprintf(fid,'\n');
    for ii = 1:NN
        fprintf(fid, '%s\t%d', these_prots{ii}, Iclusters(ii));
        for jj = 1:size(chroms,2)
            fprintf(fid,'\t%6.4f', chroms(ii,jj));
        end
        fprintf(fid,'\n');
    end
    fclose all;
    
end



%% Illustration: cluster 1. gaussians and 2. random noise

pars.best_p = 10^1;
pars.best_dens = 0;
pars.best_prec = 0.5;
pars.best_I = 10;

% 1. gaussians
Nd = 100;
Nc = 4;
sd = 0.05;
mu = [.2 .2; .8 .25; .22 .75; .68 .81];
xy = zeros(Nd,2);
for ii = 1:Nc
    xy((ii-1)*25+1 : ii*25, 1) = randn(1,25)*sd + mu(ii,1);
    xy((ii-1)*25+1 : ii*25, 2) = randn(1,25)*sd + mu(ii,2);
end
intMatrix = (1 - squareform(pdist(xy,'euclidean'))) < 0.3; 
clusters1 = clustone_mcl(intMatrix, pars, java_path);
Iclust1 = zeros(size(xy,1),1);
for ii = 1:length(clusters1)
    Iclust1(clusters1{ii}) = ii;
end

figure,hold on
scatter(xy(:,1), xy(:,2), 20)
cols = {'r' 'g' 'b' 'y' 'k' [.6 .6 .6], rand(1,3), rand(1,3), rand(1,3), rand(1,3), ...
    rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3)};
for ii = 1:length(clusters1)
    scatter(xy(clusters1{ii},1), xy(clusters1{ii},2), 25, cols{ii}, 'filled')
end
axis([0 1 0 1])


% 2. random noise
xy2 = rand(100,2)*.9 + .05;
intMatrix = (1 - squareform(pdist(xy2,'euclidean'))) < 0.3; 
clusters2 = clustone_mcl(intMatrix, pars, java_path);
Iclust2 = zeros(size(xy2,1),1);
for ii = 1:length(clusters2)
    Iclust2(clusters2{ii}) = ii;
end

figure,hold on
scatter(xy2(:,1), xy2(:,2), 20)
for ii = 1:length(clusters2)
    scatter(xy2(clusters2{ii},1), xy2(clusters2{ii},2), 25, cols{ii}, 'filled')
end
axis([0 1 0 1])


% write for R
fn = '../data/random_clusters_ex.txt';
fid = fopen(fn, 'w');
fprintf(fid, 'x\ty\tcluster\tgroup\n');
for ii = 1:length(xy)
    fprintf(fid, '%6.4f\t%6.4f\t%d\t%d\n',xy(ii,1), xy(ii,2), Iclust1(ii), 1);
end
for ii = 1:length(xy2)
    fprintf(fid, '%6.4f\t%6.4f\t%d\t%d\n',xy2(ii,1), xy2(ii,2), Iclust2(ii), 2);
end
fclose all;




