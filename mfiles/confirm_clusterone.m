% compare my code to clusterone.java for a small connection matrix MM

% connection matrices
load('/Users/Mercy/Academics/Foster/ClusterExplore/data/MM.mat')

% clusterone.java
CCj = cell(1,3);
fns = dir('E:/Greg/ClusterReliable/data/clustervalidation_*txt');
for ii = 1:length(fns)
    fid = fopen([fns(ii).folder '/' fns(ii).name]);
    cc = 0;
    while not(feof(fid))
        cc = cc+1;
        t1 = strsplit(fgetl(fid),'\t');
        CCj{ii}{cc} = nan(size(t1));
        for jj = 1:length(CCj{ii}{cc})
            CCj{ii}{cc}(jj) = str2double(t1{jj});
        end
    end
    fclose(fid);
end

% my code
CC = cell(1,length(MM));
for ii = 1:length(MM)
    CC{ii} = myclusterone(MM{ii}, 500, 0.3);
end

% visualize
figure
for ii = 1:3
    subplot(3,3,ii)
    imagesc(MM{ii})
    
    % clusterone java
    subplot(3,3,ii+3)
    Mtmp = zeros(size(MM{ii}));
    for jj = 1:length(CCj{ii})
        prots = CCj{ii}{jj};
        for kk = 1:length(prots)
            for mm = 1:length(prots)
                Mtmp(prots(kk), prots(mm)) = jj;
            end
        end
    end
    imagesc(Mtmp)
    
    % my clusterone
    subplot(3,3,ii + 6)
    Mtmp = zeros(size(MM{ii}));
    for jj = 1:length(CC{ii})
        prots = CC{ii}{jj};
        for kk = 1:length(prots)
            for mm = 1:length(prots)
                Mtmp(prots(kk), prots(mm)) = jj;
            end
        end
    end
    imagesc(Mtmp)
end
