% analyze 'data' made by ClusterExplore
% Do novel/corum complexes degrade at different rates wrt noise?


% assign a corum-like score to each non-noise complex
% just the fraction of edges that are in corum

MM = nan(10^5,2);
counter = 0;
for ii = 1:8
    ii
    %ii = 4; % dataset
    jj = 1; % noise range
    
    C0 = data.netnoise.shuffle.mcl.cluster{ii,jj};
    prots = data.Protein{ii,jj}.NoIsoform;
    
    % filter C0 to N>=5 complexes
    NN = nan(length(C0), 1);
    for mm = 1:length(C0)
        NN(mm) = length(C0{mm});
    end
    C0 = C0(NN>=5);
    NN = NN(NN>=5);
    
    corum_score = nan(length(C0),1);
    for kk = 1:length(C0)
        this_complex = prots(C0{kk});
        Ia = ismember(data.corum.proteins, this_complex);
        corum_connection_matrix = data.corum.network(Ia,Ia);
        nn_edges_in_corum = nansum(corum_connection_matrix(:)) / 2;
        
        nn = length(this_complex);
        nn_edges = nn * (nn-1) / 2;
        
        corum_score(kk) = nn_edges_in_corum / nn_edges;
    end
    
    
    % compare each non-noise complex to noised complexes
    
    % CS = maximum overlap score between original (un-noised) complex and a
    %      set of noised complexes
    CS = nan(length(corum_score), length(data.shufflenoiseRange));
    for cc = 7:8%2:length(data.shufflenoiseRange)
        C1 = data.netnoise.shuffle.mcl.cluster{ii,cc};
        % filter C1 to NN>=3
        NN = nan(length(C1), 1);
        for mm = 1:length(C1)
            NN(mm) = length(C1{mm});
        end
        C1 = C1(NN>=5);
        if length(C1)<2; continue; end
        
        % compare each original complex to each noised complex
        N0 = length(C0);
        N1 = length(C1);
        
        overlapMatrix = zeros(N0,N1);
        for mm = 1:N0
            for kk = 1:N1
                overlap = length(intersect(C0{mm},C1{kk}))^2;
                overlapMatrix(mm,kk) = overlap / length(C0{mm}) / length(C1{kk});
            end
        end
        
        CS(:,cc) = max(overlapMatrix,[],2);
    end
    
    xx = repmat(corum_score, 1, 9);
    xx = xx(:);
    yy = CS(:);
    
    I = counter+1 : counter+length(xx);
    MM(I,:) = [xx yy];
    
    counter = counter+length(xx);
end

I = not(isnan(MM(:,1))) & not(isnan(MM(:,2)));
MM = MM(I, :);

figure, hold on
scatter(MM(:,1), MM(:,2), 'filled', 'markerfacealpha', 0.1)
[R, p] = corr(MM(:,1), MM(:,2), 'type', 'spearman');
axis([0 1 0 1])
text(.5, .2, ['R=' num2str(R)])
text(.5, .1, ['R=' num2str(p)])
xlabel('Fraction of edges in CORUM, noise=0')
ylabel('Matching score, noise=0.5')
set(gcf,'units','normalized','position',[.05 .1 .5 .5],...
    'paperunits','normalized','paperposition',[.05 .1 .6 .35])
fn = 'E:\Greg\ClusterReliable\figures/novel_complexes.png';
saveas(gcf,fn)
