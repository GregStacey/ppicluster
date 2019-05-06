# Copy Figure 2 from https://arxiv.org/pdf/1706.06136.pdf

# 1000 proteins in 100 clusters
# Three intuition tests
#   1. Decays as classes are shuffled
#   2. Insensitive to number of clusters
#   3. Handles moonlighting proteins


install.packages("NMI", repos='http://cran.us.r-project.org')
install.packages("sp", repos='http://cran.us.r-project.org')
install.packages("maps", repos='http://cran.us.r-project.org')
install.packages("fossil", repos='http://cran.us.r-project.org')
install.packages("infotheo", repos='http://cran.us.r-project.org')
install.packages("FlowSOM", repos='http://cran.us.r-project.org')

source("functions.R")

fn = "/Users/gregstacey/Academics/Foster/ClusterExplore/data/suppFig01.Rda"
if (T) {
  load(fn)
} else {
  # 1. Shuffling
  tmp = seq(from=1, to=1000, by=5)
  df = data.frame(n.shuffle = tmp,
                  GA = numeric(length(tmp)),
                  MMR = numeric(length(tmp)),
                  J = numeric(length(tmp)), 
                  NMI = numeric(length(tmp)),
                  ARI = numeric(length(tmp)),
                  `F-measure` = numeric(length(tmp)),
                  MIz = numeric(length(tmp)),stringsAsFactors = F)
  clust0 = data.frame(prots = 1:1000, clusters = rep(1:100, 10))
  tmp0 = character(100)
  for (ii in 1:length(tmp0)) {
    tmp0[ii] = paste(as.character(clust0$prots[clust0$clusters==ii]), collapse=";")
  }
  for (ii in 1:nrow(df)) {
    print(round(ii/nrow(df)*100)/100)
    clust1 = clust0
    I = sample(nrow(clust0), df$n.shuffle[ii])
    clust1$clusters[I] = ceiling(runif(length(I)) * max(clust1$clusters))
    
    # convert to geomacc format
    tmp1 = character(length(unique(clust1$clusters)))
    for (jj in 1:length(tmp1)) {
      tmp1[jj] = paste(as.character(clust1$prots[clust1$clusters==jj]), collapse=";")
    }
    
    df$GA[ii] = geomacc(tmp0, tmp1)
    df$MMR[ii] = matchingratio(tmp0, tmp1)
    tmp = numeric(length(tmp0))
    for (jj in 1:length(tmp)) {
      tmp[jj] = calcA(tmp1[jj], tmp0)
    }
    df$J[ii] = mean(tmp, na.rm=T)
    df$NMI[ii] = NMI(clust0, clust1)[[1]]
    df$ARI[ii] = adj.rand.index(clust0$clusters, clust1$clusters)
    df$F.measure[ii] = FMeasure(clust0$clusters, clust1$clusters, silent=T)
    df$MIz[ii] = calcMIz(clust0$clusters,clust1$clusters,100)
  }
  df.m = melt(df, id.vars = "n.shuffle")
  df2save = list(df)
  save(df2save, file="../data/suppFig01_v02.Rda")
  
  
  # 2. Number of clusters
  tmp = seq(from=10, to=1000, by=5)
  df2 = data.frame(n.clusters = tmp,
                   GA = numeric(length(tmp)),
                   MMR = numeric(length(tmp)),
                   J = numeric(length(tmp)), 
                   NMI = numeric(length(tmp)),
                   ARI = numeric(length(tmp)),
                   `F-measure` = numeric(length(tmp)), 
                   MIz = numeric(length(tmp)), stringsAsFactors = F)
  for (ii in 1:nrow(df2)) {
    print(round(ii/nrow(df2)*100)/100)
    
    clust0 = data.frame(prots = 1:1000, clusters = rep(1:df2$n.clusters[ii], length.out=1000))
    tmp0 = character(df2$n.clusters[ii])
    for (jj in 1:length(tmp0)) {
      tmp0[jj] = paste(as.character(clust0$prots[clust0$clusters==jj]), collapse=";")
    }
    tmp0 = tmp0[!tmp0==""]
    
    clust1 = data.frame(prots = 1:1000, clusters = rep(1:df2$n.clusters[ii], length.out=1000))
    I = sample(nrow(clust0), 500)
    clust1$clusters[I] = ceiling(runif(length(I)) * max(clust1$clusters))
    # convert to geomacc format
    tmp1 = character(length(unique(clust1$clusters)))
    for (jj in 1:length(tmp1)) {
      tmp1[jj] = paste(as.character(clust1$prots[clust1$clusters==jj]), collapse=";")
    }
    tmp1 = tmp1[!tmp1==""]
    
    df2$GA[ii] = geomacc(tmp0, tmp1)
    df2$MMR[ii] = matchingratio(tmp0, tmp1)
    tmp = numeric(length(tmp0))
    for (jj in 1:length(tmp)) {
      tmp[jj] = calcA(tmp1[jj], tmp0)
    }
    df2$J[ii] = mean(tmp, na.rm=T)
    df2$NMI[ii] = NMI(clust0, clust1)[[1]]
    df2$ARI[ii] = adj.rand.index(clust0$clusters, clust1$clusters)
    df2$F.measure[ii] = FMeasure(clust0$clusters, clust1$clusters, silent=T)
    df2$MIz[ii] = calcMIz(clust0$clusters,clust1$clusters,100)
  }
  df2.m = melt(df2, id.vars = "n.clusters")
  df2save = list(df, df2)
  save(df2save, file="../data/suppFig01_v02.Rda")
  
  
  # 3. Number of moonlighting proteins
  tmp = seq(from=1, to=1000, by=5)
  df3 = data.frame(n.moonlight = tmp,
                   GA = numeric(length(tmp)),
                   MMR = numeric(length(tmp)),
                   J = numeric(length(tmp)), 
                   NMI = numeric(length(tmp)),
                   ARI = numeric(length(tmp)),
                   MIz = numeric(length(tmp)),
                   `F-measure` = numeric(length(tmp)), stringsAsFactors = F)
  prots.orig = 1:1000
  for (ii in 1:nrow(df3)) {
    print(round(ii/nrow(df3)*100)/100)
    
    clust0 = data.frame(prots = prots.orig, clusters = rep(1:100, length.out=1000))
    # prots to swap
    prots.before = sample(prots.orig, df3$n.moonlight[ii])
    I.before = clust0$prots %in% prots.before
    
    # swap them
    prots.swapped = sample(prots.orig[!prots.orig %in% prots.before], size=sum(I.before), replace = T)
    clust0$prots[I.before] = prots.swapped
    
    tmp0 = character(100)
    for (jj in 1:length(tmp0)) {
      tmp0[jj] = paste(as.character(clust0$prots[clust0$clusters==jj]), collapse=";")
    }
    
    clust1 = clust0
    tmp1 = tmp0
    
    df3$GA[ii] = geomacc(tmp0, tmp1)
    df3$MMR[ii] = matchingratio(tmp0, tmp1)
    tmp = numeric(length(tmp0))
    for (jj in 1:length(tmp)) {
      tmp[jj] = calcA(tmp1[jj], tmp0)
    }
    df3$J[ii] = mean(tmp, na.rm=T)
    df3$NMI[ii] = NMI(clust0, clust1)[[1]]
    df3$ARI[ii] = adj.rand.index(clust0$clusters, clust1$clusters)
    df3$F.measure[ii] = FMeasure(clust0$clusters, clust1$clusters, silent=T)
    df3$MIz[ii] = calcMIz(clust0$clusters,clust1$clusters,100)
  }
  df3.m = melt(df3, id.vars = "n.moonlight")
  df2save = list(df, df2, df3)
  save(df2save, file="../data/suppFig01_v02.Rda")
  
  
  # 4. Novel clusters in set 2
  tmp = seq(from=0, to=995, by=5)
  df4 = data.frame(n.novel = tmp,
                   GA = numeric(length(tmp)),
                   MMR = numeric(length(tmp)),
                   J = numeric(length(tmp)), 
                   NMI = numeric(length(tmp)),
                   ARI = numeric(length(tmp)),
                   MIz = numeric(length(tmp)),
                   `F-measure` = numeric(length(tmp)), stringsAsFactors = F)
  # clust0 as usual
  clust0 = data.frame(prots = 1:1000, clusters = rep(1:100, 10))
  tmp0 = character(100)
  for (jj in 1:length(tmp0)) {
    tmp0[jj] = paste(as.character(clust0$prots[clust0$clusters==jj]), collapse=";")
  }
  for (ii in 1:nrow(df4)) {
    print(round(ii/nrow(df4)*100)/100)
    
    # in clust1, remove some proteins
    clust1 = clust0
    I = sample(nrow(clust0), df4$n.novel[ii])
    clust1$clusters[I] = NA
    unq.clusts = unique(clust1$clusters[!is.na(clust1$clusters)])
    tmp1 = character(length(unq.clusts))
    for (jj in 1:length(tmp1)) {
      I2 = clust1$clusters==jj & !is.na(clust1$clusters)
      tmp1[jj] = paste(as.character(clust1$prots[I2]), collapse=";")
    }

    df4$GA[ii] = geomacc(tmp0, tmp1)
    df4$MMR[ii] = matchingratio(tmp1, tmp0)
    tmp = numeric(length(tmp0))
    for (jj in 1:length(tmp)) {
      tmp[jj] = calcA(tmp1[jj], tmp0)
    }
    df4$J[ii] = mean(tmp, na.rm=T)
    df4$NMI[ii] = NMI(clust0, clust1)[[1]]
    df4$ARI[ii] = adj.rand.index(clust0$clusters, clust1$clusters)
    #df4$F.measure[ii] = FMeasure(clust0.fm$clusters, clust1$clusters, silent=T)
    df4$MIz[ii] = calcMIz(clust0$clusters,clust1$clusters,100)
  }
  df4$F.measure = NA
  df4.m = melt(df4, id.vars = "n.novel")
  
  df2save = list(df, df2, df3, df4)
  save(df2save, file="../data/suppFig01_v02.Rda")
}


# 
# # A - random
# x = seq(from=0, to=1, by=0.01)
# df.fake1 = data.frame(x=x, 
#                       y=exp(-x * 3)*.9 + .1)
# ggplot(df.fake1, aes(x=x, y=y)) + geom_line() +
#   ylab("Similarity") + xlab("Fraction random in 2") + theme_bw()+
#   coord_cartesian(xlim=c(-0.02, 1.02), ylim = c(0,1.02))
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1A_v01.pdf"
# ggsave(fn,width=3, height=2.4)
# 
# df.m = melt(df, id.vars = "n.shuffle")
# ggplot(df.m, aes(x=n.shuffle/1000, y=value, color=variable)) + geom_line() +
#   ylab("Similarity") + xlab("Fraction random in 2") + theme_bw() +
#   theme(legend.position = "none") +
#   coord_cartesian(xlim=c(-0.02, 1.02), ylim = c(0,1.02))
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1A2_v01.pdf"
# ggsave(fn,width=3, height=2.4)
# 
# 
# 
# # B - number
# x = seq(from=0, to=1000, by=10)
# df.fake1 = data.frame(x=x, 
#                       y=rep(0.5, length=length(x)))
# ggplot(df.fake1, aes(x=x, y=y)) + geom_line() +
#   ylab("Similarity") + xlab("Number of clusters in 1 and 2") + theme_bw() +
#   coord_cartesian(xlim=c(-20, 1020), ylim = c(0,1.02))
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1B_v01.pdf"
# ggsave(fn,width=3, height=2.4)
# 
# df.m = melt(df2, id.vars = "n.clusters")
# ggplot(df.m, aes(x=n.clusters, y=value, color=variable)) + geom_line() +
#   ylab("Similarity") + xlab("Number of clusters in 1 and 2") + theme_bw() +
#   theme(legend.position = "none") +
#   coord_cartesian(xlim=c(-20, 1020), ylim = c(0,1.02))
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1B2_v01.pdf"
# ggsave(fn,width=3, height=2.4)
# 
# 
# 
# # C - moonlighting
# x = seq(from=0, to=1, by=0.01)
# df.fake1 = data.frame(x=x, 
#                       y=rep(1, length=length(x)))
# ggplot(df.fake1, aes(x=x, y=y)) + geom_line() +
#   ylab("Similarity") + xlab("Fraction moonlighting in 1 and 2") + theme_bw() +
#   coord_cartesian(xlim=c(-0.02, 1.02), ylim = c(0,1.02))
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1C_v01.pdf"
# ggsave(fn,width=3, height=2.4)
# 
# df.m = melt(df3, id.vars = "n.moonlight")
# ggplot(df.m, aes(x=n.moonlight/1000, y=value, color=variable)) + geom_line() +
#   ylab("Similarity") + xlab("Fraction moonlighting in 1 and 2") + theme_bw() +
#   scale_colour_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF","#00BA38")) +
#   theme(legend.position = "none") +
#   coord_cartesian(xlim=c(-0.02, 1.02), ylim = c(0,1.02))
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1C2_v01.pdf"
# ggsave(fn,width=3, height=2.4)
# 
# 
# 
# # D - novel members
# x = seq(from=0, to=1, by=0.01)
# y = exp(-x * 2)
# df.fake1 = data.frame(x=x, 
#                       y = (y - min(y)) / (max(y) - min(y)))
# ggplot(df.fake1, aes(x=x, y=y)) + geom_line() +
#   ylab("Similarity") + xlab("Fraction removed in 2") + theme_bw()+
#   coord_cartesian(xlim=c(-0.02, 1.02), ylim = c(0,1.02))
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1D_v01.pdf"
# ggsave(fn,width=3, height=2.4)
# 
# df.m = melt(df4, id.vars = "n.novel")
# ggplot(df.m, aes(x=n.novel/1000, y=value, color=variable)) + geom_line() +
#   ylab("Similarity") + xlab("Fraction removed in 2") + theme_bw() +
#   theme(legend.position = "none") +
#   coord_cartesian(xlim=c(-0.02, 1.02), ylim = c(0,1.02))
# fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1D2_v01.pdf"
# ggsave(fn,width=3, height=2.4)
