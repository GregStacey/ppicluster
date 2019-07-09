% compare my code to clusterone.java for a small connection matrix MM

% connection matrices
load('E:/Greg/ClusterReliable/data/MM.mat')

% clusterone.java
CCj = cell(1,3);
java_path = 'E:/Greg/ClusterReliable/java/cluster_one-1.0.jar';
for ii = 1:3
    CCj{ii} = clusterone_java(MM{ii}, 2, 0.3, java_path);
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


