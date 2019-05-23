
require(ggplot2)
source("functions.R")

#fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/data/fig2B_v05.Rda"
#load(fn) # sim, Ji
fn = "../data/clusters_full_netw_walktrap.txt"
Ji = as.data.frame(read_tsv(fn))
Ji = Ji[!Ji$algorithm=="hierarchical",]
Ji$measure = Ji$algorithm
Ji$algorithm[Ji$algorithm=="co_mcl"] = "CO+MCL"
Ji$algorithm[Ji$algorithm=="co"] = "CO"
Ji$algorithm[Ji$algorithm=="mcl"] = "MCL"
Ji$algorithm[Ji$algorithm=="pam"] = "k-Med"
Ji$algorithm[Ji$algorithm=="walk"] = "walktrap"

# make cluster.size, remove all clusters with size<3
Ji$cluster.size = sapply((sapply(Ji$cluster, strsplit, ";")), length)
Ji = Ji[Ji$cluster.size>2,]

# make size.factor
Ji$noise_mag = as.numeric(Ji$noise_mag)
Ji$size.factor = character(nrow(Ji))
Ji$size.factor[Ji$cluster.size<=3] = "N<=3"
Ji$size.factor[Ji$cluster.size<=6 & Ji$cluster.size>3] = "3<N<=6"
Ji$size.factor[Ji$cluster.size<=12 & Ji$cluster.size>6] = "6<N<=12"
Ji$size.factor[Ji$cluster.size>12] = "N>12"
Ji$size.factor = factor(Ji$size.factor, levels= c("N<=3", "3<N<=6", "6<N<=12", "N>12"))

# make size.prctile
Ji$size.prctile = numeric(nrow(Ji))
for (ii in 1:length(unique(Ji$noise_mag))) {
  for (jj in 1:length(unique(Ji$algorithm))) {
    for (kk in 1:length(unique(Ji$iter))) {
      I = Ji$noise_mag==sort(unique(Ji$noise_mag))[ii] & 
        Ji$algorithm==sort(unique(Ji$algorithm))[jj] & 
        Ji$iter==sort(unique(Ji$iter))[kk]
      Ji$size.prctile[I] = rank(Ji$cluster.size[I]) / sum(I)
      #Ji$size.prctile[I] = quantile(Ji$cluster.size[I], seq(from=0, to=1, length=sum(I)))
    }
  }
}

# make noise.factor
Ji$noise.factor = factor("fpr=0", levels = c("fpr=0", "fpr<=0.05","0.05<fpr<=0.20","0.2<fpr<=0.5","0.50<fpr"))
Ji$noise.factor[Ji$noise_mag<=0.05] = as.factor("fpr<=0.05")
Ji$noise.factor[Ji$noise_mag==0] = as.factor("fpr=0")
Ji$noise.factor[Ji$noise_mag>0.05 & Ji$noise_mag<=0.2] = "0.05<fpr<=0.20"
Ji$noise.factor[Ji$noise_mag>0.2 & Ji$noise_mag<=0.5] = "0.2<fpr<=0.5"
Ji$noise.factor[Ji$noise_mag>=0.5] = "0.50<fpr"

fn = "../data/clusters_Ai_vs_fdr_df_z.Rda"
load(fn) # df.z


# get averages for figures
nn = 10^3
dm = data.frame(noise_mag = numeric(nn), 
                Ji1 = numeric(nn),
                Ji2 = numeric(nn),
                algorithm = character(nn),
                size.factor = character(nn),
                stringsAsFactors = F)
cc = 0
for (ii in 1:length(unique(Ji$noise_mag))) {
  this.noise = sort(unique(Ji$noise_mag))[ii]
  for (jj in 1:length(unique(Ji$algorithm))) {
    this.algorithm = sort(unique(Ji$algorithm))[jj]
    for (kk in 1:length(unique(Ji$size.factor))) {
      this.size = sort(unique(Ji$size.factor))[kk]
      I = Ji$noise_mag==this.noise & 
        Ji$algorithm==this.algorithm & 
        Ji$size.factor==this.size
      if (sum(I)<10) next
      cc = cc+1
      dm[cc,] = c(this.noise, mean(Ji$Ji1[I], na.rm=T), mean(Ji$Ji2[I], na.rm=T), this.algorithm, this.size)
    }
  }
}
dm = dm[1:cc,]
dm = dm[!is.na(dm$size.factor),]
dm$noise_mag = as.numeric(dm$noise_mag)
dm$Ji1 = as.numeric(dm$Ji1)
dm$Ji2 = as.numeric(dm$Ji2)
tmp = sort(unique(Ji$size.factor))
dm$size.factor = tmp[as.numeric(dm$size.factor)]
dm = dm[!dm$algorithm=="hierarchical",]


Ar_null_n3 = mean(df.z$Ar[df.z$size==3],na.rm=T)
Ar_null_n6 = mean(df.z$Ar[df.z$size==6],na.rm=T)
Ar_null_n12 = mean(df.z$Ar[df.z$size==12],na.rm=T)


# 4B. Cluster size vs Ji
ggplot(Ji[Ji$noise_mag>0 & Ji$iter>1,], aes(x=log10(cluster.size), y=Ji2, color=noise.factor)) + 
  geom_point(alpha=0.05) + facet_grid(~algorithm) +
  geom_smooth(method="lm") + 
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) + 
  scale_x_continuous(breaks=log10(c(3,10,25,50,100,200)), labels = c(3,10,25,50,100,200)) +
  coord_cartesian(ylim = c(-.001,1.001), xlim=c(log10(2.6),log10(200))) + theme_bw() + 
  geom_hline(yintercept=Ar_null_n3, linetype="dashed", alpha=.65, colour="#F8766D") +
  geom_hline(yintercept=Ar_null_n6, linetype="dashed", alpha=.65, colour="#7CAE00") +
  geom_hline(yintercept=Ar_null_n12, linetype="dashed", alpha=.65, colour="#00BFC4") +
  scale_color_brewer(palette = "Spectral") +
  ylab("Similarity to iter=1 (Ji)") + 
  xlab("Cluster size")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4C_v03.png"
ggsave(fn,width=10, height=3)


# 4D. Consensus adjacency matrix examples
load("/Users/gregstacey/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr.Rda")
I0 = clusters$algorithm=="mcl" & clusters$noise_mag==0.15
clusters = clusters[I0,]
I1 = which(clusters$iter==1)

clust.good = clusters$cluster[166]
clust.bad =  clusters$cluster[78]

df.adj.good = consensus.adjmat(clust.good, clusters)
df.adj.bad = consensus.adjmat(clust.bad, clusters)
df.adj.good$iter.name = factor(paste("iter=", df.adj.good$iter, sep=""))
df.adj.bad$iter.name = factor(paste("iter=", df.adj.bad$iter, sep=""))
df.adj.good$iter.name = factor(df.adj.good$iter.name, levels=paste("iter=", 1:11, sep=""))
df.adj.bad$iter.name = factor(df.adj.bad$iter.name, levels=paste("iter=", 1:11, sep=""))

df.adj.good = df.adj.good[df.adj.good$iter %in% c(1,3,5,7,10),]
ggplot(df.adj.good, aes(prots, variable)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  facet_grid(~iter.name) + theme(legend.position="none") +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4D_v03.png"
ggsave(file=fn, width=10 * .8, height=2.3 * .8)

df.adj.bad = df.adj.bad[df.adj.bad$iter %in% c(1,3,5,7,10),]
ggplot(df.adj.bad, aes(prots, variable)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  facet_grid(~iter.name) + theme(legend.position="none") +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4E_v03.png"
ggsave(file=fn, width=10 * .8, height=2.3 * .8)

