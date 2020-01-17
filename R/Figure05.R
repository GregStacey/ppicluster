
require(ggplot2)
source("functions.R")
load("../data/broken.Rda") # broken
load("../data/test_clust_perturb_4algs.Rda") # clusters.co
co.alg = function(x) clusteroneR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar")

fn = "../data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))
corum = corum[corum$Organism=="Human",]

df = data.frame(clust = clusts,
                repJ = clusters.co[[3]]$reproducibility.J,
                size = numeric(length(clusts)),
                density = numeric(length(clusts)),
                n.outside = numeric(length(clusts)),
                seg.repJ = numeric(length(clusts)),
                seg0.repJ = numeric(length(clusts)), 
                corum.J = numeric(length(clusts)), stringsAsFactors = F)
for (ii in 1:length(clusts)) {
  print(ii)
  prots = unlist(strsplit(df$clust[ii], ";"))
  nn = length(prots)
  i.within = (ints.corum$protA %in% prots & ints.corum$protB %in% prots)
  i.without = (ints.corum$protA %in% prots | ints.corum$protB %in% prots) & !i.within
  
  df$size[ii] = nn
  df$density[ii] = sum(i.within) / (nn) / (nn-1) * 2
  df$n.outside[ii] = sum(i.without)
  
  df$seg.repJ[ii] = mean(broken$Ji[broken$id.clust==ii])# & broken$frac.scrambled<=0.99])
  df$seg0.repJ[ii] = mean(broken$Ji[broken$id.clust==ii & broken$frac.scrambled==0])
  
  df$corum.J[ii] = calcA(df$clust[ii], corum$`subunits(UniProt IDs)`)
}
df[df$density==0,] = NA


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

# seg.repJ analysis
# this is too stochastic to work :(
# things don't look consistent!
# could try just upping the number of iterations...
# NO!! seg.repJ is messed up somehow, OR it's way too strict.

# fig 5B
# breaking clusters are low density
df$good.or.bad = NA
df$good.or.bad[df$repJ>0.9 & df$size>3] = 1
df$good.or.bad[df$repJ<0.5 & df$size>3] = 0
# density
ggplot(df[!is.na(df$good.or.bad),], aes(x=factor(good.or.bad), y=density)) + 
  geom_jitter(alpha=.5, width=.25) + geom_boxplot(alpha=0) + theme_bw() +
  scale_x_discrete(labels=c("0" = "Less reproducible\n(Ji<0.5)", 
                            "1" = "More reproducible\n(Ji>0.9)")) + xlab("") +
  theme(axis.text.x = element_text(angle=35, hjust=1))
ggsave("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_5c_v01.pdf",
       width=2.5, height=3)

# stats supporting density + n.outside lead to worse clusters
summary(lm(seg.repJ ~ density + n.outside + size, data=df))
summary(lm(repJ ~ density + n.outside + size, data=df[df$size>3, ]))




