
require(scatterpie)
source("functions.R")


fnsave = "../data/clusters_wshuffle_coexp.Rda"
load(fnsave)
clusters$noise_mag = as.numeric(clusters$noise_mag)

fnsave = "../data/cluster3.txt"
clusters2 = as.data.frame(read_tsv(fnsave))
clusters2$noise_mag = as.numeric(clusters2$noise_mag)
clusters2 = clusters2[clusters2$data_type=="corum",]
clusters = rbind(clusters[,1:5], clusters2)



# 2A ##### ------------------------------------------- #####
# example of noise-adding where
#   i) clusters are obviously changing
#   ii)  measures are different

I = clusters$data_type%in%"corum" & clusters$algorithm%in%"mcl" & clusters$noise_type%in%"network_shuffle"
data = clusters[I,]
set0 = data$cluster[data$noise_mag==0]

unqnoise = unique(data$noise_mag)
for (uu in c(8)) {
  set1 = data$cluster[data$noise_mag==unqnoise[uu]]
  pci = matrix(nrow=length(set1), ncol=3)#length(set0))
  for (ii in 1:length(set1)) {
    
    # find best-matching cluster
    cluster1 = unlist(strsplit(set1[ii], ";"))
    JJ = numeric(length(set0))
    for (jj in 1:length(set0)) {
      cluster0 = unlist(strsplit(set0[jj], ";"))
      nn.jacc = length(unique(c(cluster0, cluster1)))
      JJ[jj] = length(intersect(cluster0, cluster1)) / nn.jacc 
    }
    I.max = which.max(JJ)[1]
    cluster0 = unlist(strsplit(set0[I.max], ";"))
    nn.jacc = length(unique(c(cluster0, cluster1)))
    nn.intersect = length(intersect(cluster0, cluster1))
    nn.notintersect.set1 = length(cluster1) - nn.intersect
    nn.notintersect.set0 = length(cluster0) - nn.intersect
    
    pci[ii,] = c(nn.intersect, nn.notintersect.set1, nn.notintersect.set0)
  }
  pci = as.data.frame(pci)
  pci$cluster = 1:nrow(pci)
  pci$size = unlist(lapply(lapply(lapply(set1, strsplit, ";"), "[[", 1), length))
  
  # tile x and y
  # assume that cluster with sqrt(size) needs a box of width=m*sqrt(size)
  gridx = 250
  df = pci[,]
  df = df[order(df$size, decreasing = T),]
  mm = 3
  df$tile.size = round(mm * sqrt(df$size)/5)*5 * .85 + 2
  df$y = numeric(nrow(df))
  df$x = cumsum(df$tile.size) - df$tile.size[1]*.9
  df$x = (df$x %% gridx)
  y0 = 0
  I0 = 1
  for (ii in 2:nrow(df)) {
    if (abs(df$x[ii-1] - df$x[ii]) > gridx*.65) { 
      y0 = y0 + max(df$tile.size[I0:ii])/2 + 5
      I0 = ii
    }
    df$y[ii] = y0
  }
  df = df[,c(1:3,ncol(df) + seq(from=-3, to=0, by=1))]
  df[df==0] = 10^-4
  
  # plot
  ggplot() + geom_scatterpie(aes(x=x, y=y, r=sqrt(size)), data=df, cols=names(df)[1:3]) +
    #scale_fill_manual(values=c("#4dac26","#d01c8b","#dddddd")) + # (green, red, grey) = (grey, black, white)
    scale_fill_manual(values=c("grey", "white","black")) + # (green, red, grey) = (grey, black, white)
    coord_fixed() + theme_void() + theme(legend.position="none")
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_", uu, ".pdf", sep="")
  ggsave(sf, width=10, height=5)
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_", uu, ".png", sep="")
  ggsave(sf, width=10, height=5, dpi=1000)
  
  # zoom in
  ggplot() + geom_scatterpie(aes(x=x, y=y, r=sqrt(size)), data=df, cols=names(df)[1:3]) + 
    #scale_fill_manual(values=c("#4dac26","#d01c8b","#dddddd")) +
    scale_fill_manual(values=c("grey", "white","black")) + # (green, red, grey) = (grey, black, white)
    coord_fixed() + theme_void() + theme(legend.position="none") + coord_cartesian(xlim = c(80, 130))
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_zoom_", uu, ".pdf", sep="")
  ggsave(sf, width=2, height = 5)
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_zoom_", uu, ".png", sep="")
  ggsave(sf, width=2, height = 5, dpi=1000)
}


# get full network analysis MPC revisions data
tmp = get.full.analysis()
Ji = tmp[[1]]
sim = tmp[[2]]


# simple counting: 1% shuffle affects how many clusters?
df = data.frame(nn = numeric(100), ff = numeric(100), new = numeric(100))
I = Ji$network == "CORUM" & Ji$noise_mag==0.01
unqalgs = unique(Ji$algorithm)
for (ii in 1:length(unqalgs)) {
  I2 = I & Ji$algorithm==unqalgs[ii]
  df$nn[ii] = sum(Ji$Ji1[I2] < 1, na.rm=T)
  df$ff[ii] = sum(Ji$Ji1[I2] < 1, na.rm=T) / sum(I2)
  df$new[ii] = sum(Ji$Ji1[I2] == 0)
}
df = df[1:length(unqalgs), ]


# 2B ##### ------------------------------------------- #####
# line plots of Ji vs noise - CORUM
I = sim$network=="CORUM"
ggplot(sim[I,], aes(x=noise_mag, y=Ji1)) + geom_line() + 
  facet_wrap(~algorithm, nrow=2) + theme_bw() + theme(legend.position="none") + 
  geom_point(data = Ji[Ji$network=="CORUM",], alpha=0.025, color="black") +
  ylab("Similarity to original clusters (Ji)") + xlab("Network noise level")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2B_v10.png"
ggsave(fn,width=10, height=5)

# 2C ##### ------------------------------------------- #####
# line plots of Ji vs noise - MCL
alg = "MCL"
tmp = Ji[Ji$algorithm==alg,]
I = tmp$network %in% c("DrugBank","email-Eu")
tmp = rbind(tmp, tmp[I,],tmp[I,])
tmp2 = sim[sim$algorithm==alg, ]
tmp2$Ji1[tmp2$noise_mag>0.5 & tmp2$network %in% c("BIOGRID", "BioPlex")] = NA
ggplot(tmp2, aes(x=noise_mag, y=Ji1)) + geom_line() + 
  facet_wrap(~network, nrow=1) + theme_bw() + theme(legend.position="none") + 
  geom_point(data = tmp, alpha=0.025, color="black") +
  ylab("Similarity to original clusters (Ji)") + xlab("Network noise level") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75)) +
  theme(axis.text.x = element_text(angle=25))
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2C_v10.png"
ggsave(fn,width=10, height=2.5)



# get simple stats for text
Ji = get.full.analysis()[[1]]
I = Ji$noise_mag==.01 & Ji$network=="CORUM" & Ji$iter==1
unqalg = unique(Ji$algorithm)
df = data.frame(algorithm = unqalg,
                ncluster = rep(NA, length(unqalg)),
                nchanged = rep(NA, length(unqalg)), stringsAsFactors = F)
for (ii in 1:length(unqalg)) {
  ia = I & Ji$algorithm==unqalg[ii]
  df$ncluster[ii] = sum(ia)
  df$nchanged[ii] = sum(Ji$Ji1[ia] < 1, na.rm = T)
  df$J[ii] = mean(Ji$Ji1[ia])
}
df$frac.changed = df$nchanged / df$ncluster



# do non-optimal clustering
ints.corum = as.data.frame(read_tsv("../data/interactomes/corum_pairwise.txt"))
noise.range = c(0, .01)
clusts1 = clusts3 = clusts4 = clusts5 = clusts6 = clusts7 = clusts8 = clusts9 = list()
for (uu in 1:length(noise.range)) {
  ints.shuffle = shufflecorum(ints.corum, noise.range[uu])
  unqprots = unique(c(ints.shuffle[,1], ints.shuffle[,2]))
  
  # 1. hierarchical
  x = hierarch.edge.list.format(ints.shuffle)
  tmp = stats::cutree(stats::hclust(d = x, method="average"), k = 1000)
  clusts = list()
  for (ii in 1:nclust) {
    clusts[[ii]] = unqprots[tmp == ii]
  }
  clusts1[[uu]] = sapply(clusts, FUN = function(x) paste(x, collapse=";"))
  
  # 3. MCODE
  x = graph.data.frame(ints.shuffle)
  tmp = mcode(x, vwp = 0, haircut = FALSE, fluff = FALSE, fdt = 1)
  clusts = tmp[[1]] %>% lapply(., FUN = function(x) unqprots[x])
  clusts3[[uu]] = sapply(clusts, FUN = function(x) paste(x, collapse=";"))
  
  # 4. Louvain
  if (ncol(ints.shuffle)==2) ints.shuffle$weights = 1
  tmp = louvain(ints.shuffle, 1e4)
  clusts = list()
  unqclusts = unique(tmp$cluster)
  for (ii in 1:length(unqclusts)) {
    clusts[[ii]] = unqprots[tmp$cluster == unqclusts[ii]]
  }
  clusts = clusts[!tolower(unqclusts) =="singleton"]
  clusts4[[uu]] = sapply(clusts, FUN = function(x) paste(x, collapse=";"))
  
  # 5. Leiden
  adjmat = graph_from_edgelist(as.matrix(ints.shuffle[,1:2]))
  if (ncol(ints.shuffle)==3) edge.attributes(adjmat)$weight = ints.shuffle[,3]
  this.clust = leiden(adjmat, resolution_parameter = 1)
  unqprots = names(V(adjmat))
  unqclusts = sort(unique(this.clust))
  clusts = list()
  for (ii in 1:length(unqclusts)) {
    clusts[[ii]] = unqprots[this.clust == unqclusts[ii]]
  }
  clusts5[[uu]] = sapply(clusts, FUN = function(x) paste(x, collapse=";"))
  
  # 6. walktrap
  graph.object = graph_from_edgelist(as.matrix(ints.shuffle[,1:2]), directed = F)
  if (ncol(ints.shuffle)==3) edge.attributes(graph.object)$weight = ints.shuffle[,3]
  clusts = walktrap.community(graph.object, steps = 10)
  clusts6[[uu]] = sapply(clusts, FUN = function(x) paste(x, collapse=";"))
  
  # 7. pam (k-med)
  tmp = pamclust(ints.shuffle, 1000)
  clusts = sapply(tmp, strsplit, ";")
  clusts7[[uu]] = sapply(clusts, FUN = function(x) paste(x, collapse=";"))
  
  # 8. cluster one
  tmp = unlist(clusteroneR(ints.shuffle, pp=2, density_threshold = 0.3, 
                           java_path = paste(home.dir, "/java/cluster_one-1.0.jar", sep=""),
                           output_path = co.tmpdir))
  clusts = sapply(tmp, strsplit, ";")
  clusts8[[uu]] = sapply(clusts, FUN = function(x) paste(x, collapse=";"))
  
  # 9. mcl
  # G = graph.data.frame(ints.shuffle, directed=FALSE)
  # A = as_adjacency_matrix(G,type="both", names=TRUE, sparse=FALSE)
  # tmp = MCL::mcl(A, inflation = 4, addLoops = FALSE, max.iter = 100)
  # unqprots = rownames(A)
  # clusts = list()
  # unqclusts = unique(tmp$Cluster)
  # for (ii in 1:length(unqclusts)) {
  #   clusts[[ii]] = unqprots[tmp$Cluster == unqclusts[ii]]
  # }
  # clusts9[[uu]] = sapply(clusts, FUN = function(x) paste(x, collapse=";"))
}

algs = c("hierarchical","mcode","louvain","leiden","walktrap","k-med","CO", "MCL")
df2 = data.frame(algorithm = algs,
                 noise_mag = rep(.01, 8),
                 JJ = rep(NA, 8), stringsAsFactors = F)
allclusts = list(clusts1,clusts3,clusts4,clusts5,clusts6,clusts7,clusts8,clusts9)
for (ii in 1:7) {
  print(ii)
  df2$JJ[ii] = mean(calcJ2(allclusts[[ii]][[1]], allclusts[[ii]][[2]]))
}

# # test - small clusters
# I = Ji$network == "CORUM" & Ji$noise_mag==0.02 & Ji$cluster.size<=5
# agg = aggregate(Ji$Ji1[I],
#                 by = list(Ji$algorithm[I]),
#                 FUN = mean)
# print(paste(min(agg$x), " , ", max(agg$x)))
# I = Ji$network == "CORUM" & Ji$noise_mag==0.02 & Ji$cluster.size%in%(6:20)
# agg = aggregate(Ji$Ji1[I],
#                 by = list(Ji$algorithm[I]),
#                 FUN = mean)
# print(paste(min(agg$x), " , ", max(agg$x)))
# I = Ji$network == "CORUM" & Ji$noise_mag==0.02 & Ji$cluster.size%in%(21:1000)
# agg = aggregate(Ji$Ji1[I],
#                 by = list(Ji$algorithm[I]),
#                 FUN = mean)
# print(paste(min(agg$x), " , ", max(agg$x)))
# 
# I = Ji$network == "CORUM" & Ji$algorithm=="CO"
# cor.test(Ji$cluster.size[I], Ji$Ji1[I])
# 
# # 2C ##### ------------------------------------------- #####
# # violin plots of Ji1 vs noise
# 
# I = (Ji$noise_mag %in% 0 | Ji$noise_mag %in% 0.01 | Ji$noise_mag %in% 0.02) & Ji$network=="CORUM"
# ggplot(Ji[I,], aes(factor(noise_mag), y=Ji1)) + 
#   geom_violin() + geom_jitter(width=.02,alpha=.05) + facet_grid(~algorithm) +
#   ylab("Similarity to un-noised clusters (Ji)") + xlab("Network noise level") + 
#   coord_cartesian(ylim=c(0,1)) + theme_bw()
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2C_v06.pdf"
# ggsave(fn,width=10, height=2.6)
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2C_v06.png"
# ggsave(fn,width=10, height=2.6)
# 
# 
# # 2D ##### ------------------------------------------- #####
# # noise-amplification
# 
# # I = sim$measure=="ai"& sim$network=="CORUM"
# # sim$noise.amplification = (1-sim$Ji1) / sim$noise_mag
# # ggplot(sim[I,], aes(x=noise_mag, y=(noise.amplification))) + geom_line() + 
# #   facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
# #   geom_hline(yintercept = 1, linetype="dashed") +
# #   ylab("Error amplification by clustering") + xlab("Interactome FPR") +
# #   scale_y_continuous(breaks = c(1,10,20),
# #                      labels=c("1x","10x", "20x"))
# 
# sf = "../data/dfchange_02.Rda"
# load(sf)
# 
# ggplot(df.change, aes(x=n.shuffled.interactions, y=n.rearranged.edges, shape=algorithm)) + 
#   geom_point(alpha=0.6) + geom_abline(linetype="dashed") +theme_bw() +
#   xlab("Number of shuffled network edges") + 
#   ylab("Number of rearranged cluster edges") + 
#   coord_cartesian(ylim=c(0,80000)) + theme(legend.position="none")
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_3D_v01.pdf"
# ggsave(fn,width=4, height=2.5)
# 
# ggplot(df.change, aes(x=noise, y=n.rearranged.edges/n.shuffled.interactions, color=algorithm)) + 
#   geom_line() + geom_abline(linetype="dashed") +theme_bw() +
#   xlab("Network noise level") + 
#   ylab("Ratio of rearranged cluster edges\nto shuffled network edges") + theme(legend.position="none") +
#   scale_colour_grey()
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2D2_v04.pdf"
# ggsave(fn,width=4, height=3.2)
# 
# 
# # 2E,F ##### ------------------------------------------- #####
# 
# I = sim$measure=="ai" & sim$measure=="ai" & sim$network %in% c("DrugBank","email-Eu") & sim$algorithm=="MCL"
# ggplot(sim[I,], aes(x=noise_mag, y=Ji1)) + geom_line() + 
#   theme_bw() + theme(legend.position="none") + facet_grid(~network) +
#   geom_point(data = Ji[Ji$network%in% c("DrugBank","email-Eu") &Ji$algorithm=="MCL",], alpha=0.1, color="black") +
#   ylab("Similarity to un-noised clusters (Ji)") + xlab("Network noise level")
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2E_v02.png"
# ggsave(fn,width=4.2, height=2.6)
# 
# I = (Ji$noise_mag %in% 0 | Ji$noise_mag %in% 0.01 | Ji$noise_mag %in% 0.02) & Ji$network%in% c("DrugBank","email-Eu") & 
#   Ji$algorithm=="MCL"
# ggplot(Ji[I,], aes(factor(noise_mag), y=Ji1)) + 
#   geom_violin() + geom_jitter(width=.02,alpha=.1) + facet_grid(~network) +
#   ylab("Similarity to un-noised (Ji)") + xlab("Network noise level") + 
#   coord_cartesian(ylim=c(0,1)) + theme_bw()
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2F_v02.pdf"
# ggsave(fn,width=4.2, height=2.6)
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2F_v02.png"
# ggsave(fn,width=4.2, height=2.6)
# 
# 
# 
# 
# # supp of drugbank and email sets ##### ------------------------------------------- #####
# for (ss in 1:2) {
#   sets = c("DrugBank","email-Eu")
#   # line plots of J_i vs noise
#   I = sim$measure=="ai" & sim$measure=="ai" & sim$network==sets[ss]
#   ggplot(sim[I,], aes(x=noise_mag, y=Ji1)) + geom_line() + 
#     facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
#     geom_point(data = Ji[Ji$network==sets[ss],], alpha=0.025, color="black") +
#     ylab("Similarity to un-noised clusters (Ji)") + xlab("Network noise level")
#   fn = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1A",ss,"_v02.png", sep="")
#   ggsave(fn,width=10, height=2.6)
#   
#   # violin plots of Ji1 vs noise
#   I = (Ji$noise_mag %in% 0 | Ji$noise_mag %in% 0.01 | Ji$noise_mag %in% 0.02) & Ji$network==sets[ss]
#   ggplot(Ji[I,], aes(factor(noise_mag), y=Ji1)) + 
#     geom_violin() + geom_jitter(width=.02,alpha=.05) + facet_grid(~algorithm) +
#     ylab("Similarity to un-noised clusters (Ji)") + xlab("Network noise level") + 
#     coord_cartesian(ylim=c(0,1)) + theme_bw()
#   fn = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1B",ss,"_v02.png", sep="")
#   ggsave(fn,width=10, height=2.6)
#   fn = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1B",ss,"_v02.pdf", sep="")
#   ggsave(fn,width=10, height=2.6)
#   
#   # noise-amplification
#   I = sim$measure=="ai"& sim$network==sets[ss]
#   sim$noise.amplification = (1-sim$Ji1) / sim$noise_mag
#   ggplot(sim[I,], aes(x=noise_mag, y=(noise.amplification))) + geom_line() + 
#     facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
#     geom_hline(yintercept = 1, linetype="dashed") +
#     ylab("Error amplification by clustering") + xlab("Network noise level") +
#     scale_y_continuous(breaks = c(1,10,20),
#                        labels=c("1x","10x", "20x"))
#   fn = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1C",ss,"_v02.pdf", sep="")
#   ggsave(fn,width=10, height=2.6)
# }
# 
# 
# 
# # non-optimal parameters
# # ##### ------------------------------------------- #####
# fn = "../data/clusters_wshuffle_notoptimal.txt"
# Jin = as.data.frame(read_tsv(fn))
# unqalgs = unique(Ji$algorithm)
# unqmags = unique(Ji$noise_mag)
# Jin$Ji1 = numeric(nrow(Jin))
# for (ii in 1:length(unqalgs)) {
#   I0 = Ji$algorithm==unqalgs[ii] & Ji$noise_mag==0
#   ref.clusters = Jin$cluster[I0]
#   for (jj in 1:length(unqmags)) {
#     print(paste("        noise",unqmags[jj]))
#     I = which(Ji$algorithm==unqalgs[ii] & Ji$noise_mag==unqmags[jj])
#     these.clusters = Jin$cluster[I]
#     for (mm in 1:length(I)) {
#       Jin$Ji1[I[mm]] = calcA(these.clusters[mm], ref.clusters)
#     }
#   }
# }
# 
# for (ii in 1:length(unqalgs)) {
#   I = Jin$noise_mag==0.01 & Jin$algorithm==unqalgs[ii]
#   x = (mean(Jin$Ji1[I]))
#   I = Ji$noise_mag==0.01 & Ji$algorithm==unqalgs[ii] & Ji$network=="CORUM"
#   x2 = (mean(Ji$Ji1[I]))
#   print(paste(x,x2))
# }
# 
# 
# 
# 
# # small-clusters-are-bad in this noise-nonoise data
# # (not sure about this)
# # ##### ------------------------------------------- #####
# I = Ji$network=="CORUM"
# summary(lm(Ji1 ~ noise_mag + algorithm + cluster.size, data = Ji[I,]))
# 
# ggplot(Ji[I,], aes(x=log10(cluster.size), y=Ji1)) + 
#   facet_grid(noise_mag~algorithm) + 
#   geom_point(alpha=.3) + geom_smooth(method="lm")
# 







