
% load adjacency matrices
adjacency = cell(1,3);

% corum
fn = 'E:/Greg/ClusterReliable/data/allComplexes.txt';
[adjacency{1}, proteins, complexes] = corum2network(fn);
% sort by appearance in complexes
ia = nan(size(proteins));
for ii = 1:length(proteins)
    I = find(~cellfun(@isempty,strfind(complexes,proteins{ii})));
    ia(ii) = I(1);
end
[x,Isort] = sort(ia);
proteins = proteins(Isort);
adjacency{1} = adjacency{1}(Isort, Isort);

% drugbank and email
fns{1} = '';
fns{2} = 'E:\Greg\ClusterReliable\data/ChCh-Miner_durgbank-chem-chem.tsv'; % drug-drug network
fns{3} = 'E:\Greg\ClusterReliable\data/email-Eu-core.txt'; % email network
seps = {'' '\t' ' '};
for uu = 2:3
    
    % Read raw prot-chem file
    fid = fopen(fns{uu},'r');
    fgetl(fid); % header
    
    % read body of file
    cc = 0;
    edges = cell(10^6,2);
    while ~feof(fid)
        t1 = strsplit(fgetl(fid),seps{uu});
        cc = cc+1;
        edges{cc,1} = t1{1};
        edges{cc,2} = t1{2};
    end
    edges = edges(1:cc,:);
    fclose(fid);
    
    % All unique proteins
    nodes = unique(edges(1:cc, :));
    
    % Split complexes into pairwise list
    adjacency{uu} = zeros(length(nodes));
    [x,ia] = ismember(edges(:,1), nodes);
    [x,ib] = ismember(edges(:,2), nodes);
    for ii = 1:length(ia)
        adjacency{uu}(ia(ii),ib(ii)) = 1;
        adjacency{uu}(ib(ii),ia(ii)) = 1;
    end
end

% sort drugbank by optimal clustering
[~, clusts] = mymcl(adjacency{2}, 32, 4, 10^-12, 20);
Isort = nan(size(adjacency{2},1),1);
for ii = 1:length(clusts)
    Isort(clusts{ii}) = ii;
end
Isort(isnan(Isort)) = length(clusts)+1;
[x,Isort] = sort(Isort);
adjacency{2} = adjacency{2}(Isort, Isort);


% sort email by department affiliation
fn = 'E:/Greg/ClusterReliable/data/email-Eu-core-department-labels.txt';
fid = fopen(fn);
node = nan(size(adjacency{3},1),1);
department_label = nan(size(adjacency{3},1),1);
cc = 0;
while not(feof(fid))
    t1 = strsplit(fgetl(fid), ' ');
    cc = cc+1;
    node(cc) = str2double(t1{1}) + 1;
    department_label(cc) = str2double(t1{2});
end
[x,Isort] = sort(department_label);
adjacencey{3} = adjacency{3}(Isort, Isort);


for ii = 2:3
    figure
    imagesc(adjacency{ii})
    caxis([0 1])
    axis square
    colormap bone
    set(gca,'xtick',0:20:80,'xticklabel','','ytick',0:20:80,'yticklabel','');
    set(gcf,'paperunits','inches','paperposition',[1 1 8.9*.75 8.9*.75],...
        'units','inches','position',[1 1 8.9*.75 8.9*.75])
    sf = ['E:/Greg/ClusterReliable/figures/network_adjacency' num2str(ii)];
    print(sf, '-dpng', '-r1000');
end

