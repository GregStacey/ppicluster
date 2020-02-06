
# make figure showing that breakage is reproducible


require(ggplot2)
source("functions.R")
load("../data/broken.Rda") # broken
load("../data/test_clust_perturb_4algs.Rda") # clusters.co
co.alg = function(x) clusteroneR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar")

fn = "../data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))
corum = corum[corum$Organism=="Human",]



# tool figure

# 1. motivation - breakage is correlated between noise iterations
if (1){
  load("../data/clust-iters.Rda")
} else {
  clust.iters1 = clust.perturb(ints.corum, clustering.algorithm = co.alg, noise=0.01, iters=10)
  clust.iters2 = clust.perturb(ints.corum, clustering.algorithm = co.alg, noise=0.05, iters=10)
  clust.iters3 = clust.perturb(ints.corum, clustering.algorithm = co.alg, noise=0.1, iters=10)
  clust.iters4 = clust.perturb(ints.corum, clustering.algorithm = co.alg, noise=0.2, iters=10)
  clust.iters5 = clust.perturb(ints.corum, clustering.algorithm = co.alg, noise=0.5, iters=10)
  save(clust.iters1,clust.iters2,clust.iters3,clust.iters4,clust.iters5, file = "../data/clust-iters.Rda")
}
jmat = matrix(nrow=length(clusts), ncol=10)
for (ii in 1:length(clusts)) {
  jmat[ii,] = as.numeric(unlist(strsplit(clust.iters3$all.iterations.repJ[ii], ";")))
}

# fig 5A
# scatter of individual replicates
tmp = as.data.frame(jmat)
tmp$size = df$size
ggplot(tmp, aes(x=V4, y=V10, size=size)) + 
  geom_jitter(alpha=.4, height = .015, width=.015) +
  theme_bw() + xlab("Ji (iteration 1)\n10% noise") + ylab("Ji (iteration 2)\n10% noise") +
  theme(legend.position = "none")
ggsave("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_5a_v01.pdf",
       width=3, height=3)

# fig 5B
# scatter of the mean of replicates
tmp = data.frame(xx = clust.iters2$reproducibility.J, 
                 yy = clust.iters3$reproducibility.J, size=df$size)
ggplot(tmp, aes(x=yy, y=xx, size=size)) + 
  geom_point(alpha=.25) + theme_bw() + 
  scale_size(guide = "legend",breaks = c(3, 25, 50, 100)) +
  xlab("Ji (10 iterations, average)\n10% noise") + 
  ylab("Ji (10 iterations, average)\n5% noise")
ggsave("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_5b_v01.pdf",
       width=3.85, height=3)



# fig 5c
# show single replicate correlation for all algorithms
###
load("../data/enrichment_clustperturb_kmed.Rda")
tmp = clusters.kmed
load("../data/enrichment_clustperturb_mcl.Rda")
tmp2 = clusters.mcl
load("../data/enrichment_clustperturb.Rda")
clusters.kmed = tmp
clusters.mcl = tmp2
rm(tmp)
rm(tmp2)
clusts = list("clusters.kmed" = clusters.kmed,
                "clusters.mcl" = clusters.mcl,
                "clusters.walk" = clusters.walk,
                "clusters.co" = clusters.co)
algs = c("k-Med","MCL", "walktrap", "CO")
###
jmat = list()
for (ii in 1:length(clusts)) {
  clusts[[ii]]$size = sapply((sapply(clusts[[ii]]$cluster, strsplit, ";")), length)
  jmat[[ii]] = matrix(nrow=nrow(clusts[[ii]]), ncol=100)
  for (jj in 1:nrow(clusts[[ii]])) {
    jmat[[ii]][jj,] = as.numeric(unlist(strsplit(clusts[[ii]]$all.iterations.repJ[jj], ";")))
  }
}
df2 = data.frame(algorithm = character(1e6),
                rr = numeric(1e6), stringsAsFactors = F)
cc = 0
for (uu in 1:length(algs)) {
  tmp = cor(jmat[[uu]])
  tmp = tmp[lower.tri(tmp, diag = F)]
  I = (cc+1) : (cc+length(tmp))
  df2$algorithm[I] = algs[uu]
  df2$rr[I] = tmp
  cc = cc+length(I)
}
df2 = df2[1:cc,]

ggplot(df2, aes(x=algorithm, y=rr)) + 
  geom_boxplot(alpha = 0) + geom_jitter(alpha=.002, width = .25) + 
  theme_bw() + xlab("Clustering algorithm") +
  ylab("Noise iteration pairs,\nJi correlation (Pearson)")
ggsave("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_5c_v01.png",
       width=3, height=3)



# fig 5d
# breaking clusters are low density
load("../data/clusters_with_groundtruth.Rda")
cc =0
df = data.frame(repJ = numeric(1e5),
                algorithm = numeric(1e5),
                size = numeric(1e5),
                density = numeric(1e5),
                corum.J = numeric(1e5), stringsAsFactors = F)
algs = c("k-Med", "MCL", "walktrap", "CO")
for (ii in 1:length(clusters)) {
  for (jj in 1:nrow(clusters[[ii]])) {
    print(jj)
    prots = unlist(strsplit(clusters[[ii]]$cluster[jj], ";"))
    nn = length(prots)
    i.within = (ints.corum$protA %in% prots & ints.corum$protB %in% prots)
    #i.without = (ints.corum$protA %in% prots | ints.corum$protB %in% prots) & !i.within
    #df$n.outside[ii] = sum(i.without)
    
    #df$seg.repJ[ii] = mean(broken$Ji[broken$id.clust==ii])# & broken$frac.scrambled<=0.99])
    #df$seg0.repJ[ii] = mean(broken$Ji[broken$id.clust==ii & broken$frac.scrambled==0])
    
    cc = cc+1
    df$size[cc] = nn
    df$density[cc] = sum(i.within) / (nn) / (nn-1) * 2
    df$corum.J[cc] = clusters[[ii]]$corumJ[jj]
    df$repJ[cc] = clusters[[ii]]$reproducibility.J[jj]
    df$algorithm[cc] = algs[ii]
  }
}
df = df[1:cc,]
df = df[df$density>0, ]
df$good.or.bad = rep(NA, nrow(df))
df$good.or.bad[df$repJ<0.5] = "bad"
df$good.or.bad[df$repJ>0.75] = "good"

# density
ggplot(df[!is.na(df$good.or.bad),], aes(x=factor(algorithm), color=good.or.bad, fill=good.or.bad, y=density)) + 
  geom_jitter(alpha=.5, width=.25) + geom_boxplot(alpha=0) + theme_bw() +
  scale_x_discrete(labels=c("0" = "Less reproducible\n(Ji<0.5)", 
                            "1" = "More reproducible\n(Ji>0.9)")) + xlab("") +
  theme(axis.text.x = element_text(angle=35, hjust=1))
ggsave("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_5d_v01.pdf",
       width=2.5, height=3)

ggplot(df, aes(x=density, y=repJ, color=algorithm, size=size)) + geom_point(alpha=.4) +
  theme_bw() + geom_smooth(se=F, method="lm")

alg = "k-Med"
cor.test(df$density[df$algorithm==alg], df$repJ[df$algorithm==alg])

