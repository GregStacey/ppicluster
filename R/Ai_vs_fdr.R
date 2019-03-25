

# define functions
calcA = function(this.cluster, clusters) {
  # assignment reproducibility
  # essentially the maximum Jaccard index
  #
  # this.cluster = single string, with semicolon-delimited IDs
  # clusters = vector of strings, all semicolon-delimited IDs
  
  JJ = numeric(length(clusters))
  for (ii in 1:length(clusters)) {
    that.cluster = unlist(strsplit(clusters[ii], ";"))
    JJ[ii] = length(intersect(this.cluster, that.cluster)) / 
      length(unique(c(this.cluster, that.cluster)))
  }
  A = max(JJ, na.rm=T)
  return(A)
}



# get clusters
fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr.txt"
data = read_tsv(fn)
clusters = data[data$iter==1,]
clusters$size = unlist(lapply(sapply(clusters$cluster, FUN=strsplit, ";"), length))
clusters$size.factor = character(nrow(clusters))
clusters$size.factor[clusters$size<=3] = "N<=3"
clusters$size.factor[clusters$size<=6 & clusters$size>3] = "3<N<=6"
clusters$size.factor[clusters$size<=12 & clusters$size>6] = "6<N<=12"
clusters$size.factor[clusters$size>12] = "N>12"
clusters$size.factor = factor(clusters$size.factor, levels= c("N<=3", "3<N<=6", "6<N<=12", "N>12"))


# all unique proteins
rdntprots = unlist(sapply(clusters$cluster, FUN=strsplit, ";"))
unqprots = unique(rdntprots)


# calculate Ai for each cluster
I.order = sample(nrow(clusters), nrow(clusters))
clusters$Ar = numeric(nrow(clusters))
clusters$Ar.zscore = numeric(nrow(clusters))
for (ii in 1:nrow(clusters)) {
  print(ii)
  this.cluster = unlist(strsplit(clusters$cluster[I.order[ii]], ";"))
  this.algorithm = clusters$algorithm[I.order[ii]]
  this.fdr = clusters$noise_mag[I.order[ii]]
  
  # calculate Ai, complex reproducibility
  unqIters = 2:10
  Ar = numeric(length(unqIters))
  #Ar_z_iter = 2
  #Ar_null = matrix(nrow = length(unqIters), ncol = Ar_z_iter)
  for (jj in 1:length(unqIters)) {
    I = data$algorithm %in% this.algorithm & data$noise_mag==this.fdr & data$iter==unqIters[jj]
    Ar[jj] = calcA(this.cluster, data$cluster[I])
    
    # calculate null distribution
    #for (kk in 1:Ar_z_iter) {
    #  cluster0.fake = unqprots[sample(length(unqprots), length(cluster0))]
    #  Ar_null[jj,kk] = calcA(cluster0.fake, data$cluster[I])
    #}
  }

  clusters$Ar[I.order[ii]] = mean(Ar, na.rm=T)
  clusters$Ar.zscore[I.order[ii]] = (mean(Ar, na.rm=T) - mean(Ar_null, na.rm=T)) / sd(Ar_null, na.rm=T)
}


# calculate NULL Ai
sizeRange = c(3,6,15)
iterMax = 25
df.z = data.frame(size = numeric(10^4), Ar = numeric(10^4), 
                  iter = numeric(10^4), algorithm = numeric(10^4),
                  noise_mag = numeric(10^4), stringsAsFactors = F)
unqAlgs = unique(clusters$algorithm)
unqNoise = unique(clusters$noise_mag)
cc = 0
for (ii in 1:length(sizeRange)) {
  this.size = sizeRange[ii]
  for (jj in 1:length(unqAlgs)) {
    this.alg = unqAlgs[jj]
    for (kk in 1:length(unqNoise)) {
      this.noise = unqNoise[kk]
      I = clusters$algorithm%in%this.alg & clusters$noise_mag==this.noise
      for (mm in 1:iterMax) {
        cc = cc+1
        print(cc)
        cluster0.fake = rdntprots[sample(length(rdntprots), this.size)]
        Ar_null = calcA(cluster0.fake, data$cluster[I])
        df.z$size[cc] = this.size
        df.z$Ar[cc] = Ar_null
        df.z$iter[cc] = mm
        df.z$algorithm[cc] = this.alg
        df.z$noise_mag[cc] = this.noise
      }
    }
  }
}
df.z = df.z[1:cc,]
Ar_null_n3 = mean(df.z$Ar[df.z$size==3],na.rm=T)
Ar_null_n15 = mean(df.z$Ar[df.z$size==15],na.rm=T)


# get averages for figures
nn = 10^3
dm = data.frame(x = numeric(nn), 
                y = numeric(nn),
                algorithm = character(nn),
                size.factor = character(nn),
                stringsAsFactors = F)
cc = 0
for (ii in 1:length(unique(clusters$noise_mag))) {
  this.noise = sort(unique(clusters$noise_mag))[ii]
  for (jj in 1:length(unique(clusters$algorithm))) {
    this.algorithm = sort(unique(clusters$algorithm))[jj]
    for (kk in 1:length(unique(clusters$size.factor))) {
      this.size = sort(unique(clusters$size.factor))[kk]
      cc = cc+1
      I = clusters$noise_mag==this.noise & 
        clusters$algorithm%in%this.algorithm & 
        clusters$size.factor%in%this.size
      dm[cc,] = c(this.noise, mean(clusters$Ar[I], na.rm=T), this.algorithm, this.size)
    }
  }
}
dm = dm[1:cc,]
dm$x = as.numeric(dm$x)
dm$y = as.numeric(dm$y)
tmp = sort(unique(clusters$size.factor))
dm$size.factor = tmp[as.numeric(dm$size.factor)]

ggplot(clusters, aes(x=noise_mag, y=Ar, color=algorithm)) + 
  geom_point(alpha=0.25) + facet_wrap(~size.factor) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_line(data=dm, aes(x=x, y=y,color=algorithm), size=2, alpha=.6)+
  geom_hline(yintercept=Ar_null_n3, linetype="dashed", alpha=.8)+
  geom_hline(yintercept=Ar_null_n15, linetype="dashed", alpha=.8) + 
  annotate("text",x=.7,y=Ar_null_n3+.032, label="N=3") +
  annotate("text",x=.7,y=Ar_null_n15+.032, label="N=15")
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/AvsFDR_algorithm.png")

ggplot(clusters, aes(x=noise_mag, y=Ar, color=size.factor)) + 
  geom_point(alpha=0.25) +  facet_wrap(~algorithm) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") + 
  geom_line(data=dm, aes(x=x, y=y,color=size.factor), size=2, alpha=.6) +
  coord_cartesian(ylim=c(0,1)) +
  geom_hline(yintercept=Ar_null_n3, linetype="dashed", alpha=.65)+
  geom_hline(yintercept=Ar_null_n15, linetype="dashed", alpha=.65) + 
  annotate("text",x=.7,y=Ar_null_n3+.015, label="N=3", alpha=.75) +
  annotate("text",x=.7,y=Ar_null_n15+.015, label="N=15", alpha=.75)
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/AvsFDR_size.png")

# zoom in on fdr=0.01
I = clusters$noise_mag %in% 0 | clusters$noise_mag %in% 0.01 | clusters$noise_mag %in% 0.02
ggplot(clusters[I,], aes(factor(noise_mag), y=Ar)) + 
  geom_violin() + geom_jitter(width=.02,alpha=.05) + facet_grid(algorithm ~ size.factor) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") + 
  coord_cartesian(ylim=c(0,1))
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/AvsFDR_violin_fdr01.png")

# single panels for size <=3, 5-7, and >13
I = clusters$size<=3
ggplot(clusters[I,], aes(x=noise_mag, y=Ar, color=algorithm)) +
  geom_point(alpha=0.25) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") +
  coord_cartesian(ylim=c(0,1)) +
  geom_line(data=dm[dm$size.factor%in%"N<=3",], aes(x=x, y=y,color=algorithm), size=2, alpha=.6) +
  geom_hline(yintercept=mean(df.z$Ar[df.z$size==3], na.rm=T), linetype="dashed")
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/AvsFDR_N3.png")

I = clusters$size>=5 & clusters$size<=7
ggplot(clusters[I,], aes(x=noise_mag, y=Ar, color=algorithm)) +
  geom_point(alpha=0.25) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") +
  coord_cartesian(ylim=c(0,1)) +
  geom_line(data=dm[dm$size.factor%in%"3<N<=6",], aes(x=x, y=y,color=algorithm), size=2, alpha=.6) +
  geom_hline(yintercept=mean(df.z$Ar[df.z$size==6], na.rm=T), linetype="dashed")
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/AvsFDR_N6.png")

I = clusters$size>=13
ggplot(clusters[I,], aes(x=noise_mag, y=Ar, color=algorithm)) +
  geom_point(alpha=0.25) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") +
  coord_cartesian(ylim=c(0,1)) +
  geom_line(data=dm[dm$size.factor%in%"N>12",], aes(x=x, y=y,color=algorithm), size=2, alpha=.6) +
  geom_hline(yintercept=mean(df.z$Ar[df.z$size==15], na.rm=T), linetype="dashed")
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/AvsFDR_N10.png")

ggplot(df.z, aes(x=Ar, fill=factor(size))) + geom_density()
