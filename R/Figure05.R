
require(ggplot2)
require(igraph)
require(ggnetwork)
require(intergraph)
source("/Users/Mercy/Academics/Foster/ClusterExplore/R/functions.R")

load("/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr.Rda")
load("/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr_df_z.Rda")



# hairball consensus clusters
I.iter1 = clusters$iter==1
unqIters = 2:10
unqalg = unique(clusters$algorithm)
unqnoise = sort(unique(clusters$noise_mag))
ii = 1 # algorithm

for (jj in c(9,3,5,7)) {
  print(unqnoise[jj])
  I0 = clusters$algorithm %in% unqalg[ii] & clusters$noise_mag==unqnoise[jj]
  I = which(I.iter1 & I0)
  
  # pick a cluster
  for (kk in 1:length(I)) {
    print(paste("   ", kk))
    this.cluster = unlist(strsplit(clusters$cluster[I[kk]], ";"))
    this.size = length(this.cluster)
    this.size = formatC(this.size, width=3, flag="0")
    if (length(this.cluster)>150) next
    
    # find its best match in the other iters
    bestMatches = character(length(unqIters))
    for (mm in 1:length(unqIters)) {
      set0 = clusters$cluster[I0 & clusters$iter==unqIters[mm]]
      JJ = numeric(length(set0))
      for (uu in 1:length(set0)) {
        that.cluster = unlist(strsplit(set0[uu], ";"))
        JJ[uu] = length(intersect(this.cluster, that.cluster)) / 
          length(unique(c(this.cluster, that.cluster)))
      }
      I.match = which.max(JJ)
      bestMatches[mm] = set0[I.match]
    }
    
    # make consensus adjacency matrix
    allProts = unique(c(this.cluster, unlist(lapply(bestMatches, strsplit, ";"))))
    adjmat = matrix(numeric(length(allProts)^2), nrow=length(allProts), ncol=length(allProts))
    df.adjmat = data.frame(prots = character(0),
                           variable = character(0),
                           value = numeric(0), iter=numeric(0), stringsAsFactors = F)
    bestMatches = c(paste(this.cluster,collapse=";"), bestMatches)
    Ar = numeric(length(bestMatches))
    for (mm in 1:length(bestMatches)) {
      that.cluster = unlist(strsplit(bestMatches[mm], ";"))
      for (uu in 1:length(that.cluster)) {
        ia = which(allProts %in% that.cluster[uu])
        for (vv in 1:length(that.cluster)) {
          if (uu>=vv) next
          ib = which(allProts %in% that.cluster[vv])
          adjmat[ia,ib] = adjmat[ia,ib]+1
          adjmat[ib,ia] = adjmat[ib,ia]+1
        }
      }
      
      Ar[mm] = calcA(this.cluster, paste(that.cluster, collapse=";"))
      df.tmp = as.data.frame(adjmat)
      df.tmp$prots = allProts
      df.tmp2 = melt(df.tmp, id.vars="prots")
      df.tmp2$iter = mm
      df.tmp2$value = df.tmp2$value * length(bestMatches) / mm
      df.adjmat = rbind(df.adjmat, df.tmp2)
    }
    adjmat = adjmat / length(unqIters)

    # hairball
    tryCatch({
      # plot it
      p1 = hairball(adjmat, this.cluster, allProts)
      fn = paste("/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/hairballs/",
                 unqalg[ii], "_", unqnoise[jj], "_", kk, ".png", sep="")
      ax = 3 * sqrt(length(allProts) / 15)
      ggsave(p1, file=fn, width=ax, height=ax)
    },error=function(cond) {
    },warning=function(cond) {
    },finally={
    })
    
    # adjacency matrix
    if (length(allProts)>4 & length(allProts)<100) {
      df.adjmat$prots = factor(df.adjmat$prots, levels = allProts)
      p2 = ggplot(df.adjmat, aes(prots, variable)) +
        geom_tile(aes(fill = value), color = "white") +
        scale_fill_gradient(low = "white", high = "steelblue") + 
        facet_wrap(~iter) + theme(legend.position="none") +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) + 
        annotate("text", x = length(allProts)-4, y = length(allProts)-2, 
                 label = round(mean(Ar) * 1000)/1000)
      fn = paste("/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/adjacency/","size",
                 this.size, "_",unqalg[ii], "_", unqnoise[jj], "_", kk, ".png", sep="")
      ggsave(p2, file=fn, width=8, height=7)
    }
  }
}




# 4 ##### ------------------------------------------- #####
# % reproducible vs FPR

load("/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr.Rda")
load("/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr_df_z.Rda")

clusters$noise_mag = as.numeric(clusters$noise_mag)
clusters$size = unlist(lapply(sapply(clusters$cluster, FUN=strsplit, ";"), length))
clusters$size.factor = character(nrow(clusters))
clusters$size.factor[clusters$size<=3] = "N<=3"
clusters$size.factor[clusters$size<=6 & clusters$size>3] = "3<N<=6"
clusters$size.factor[clusters$size<=12 & clusters$size>6] = "6<N<=12"
clusters$size.factor[clusters$size>12] = "N>12"
clusters$size.factor = factor(clusters$size.factor, levels= c("N<=3", "3<N<=6", "6<N<=12", "N>12"))

# get averages for figures
nn = 10^3
dm = data.frame(x = numeric(nn), 
                y = numeric(nn),
                y2 = numeric(nn),
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
        clusters$size.factor%in%this.size & clusters$iter==1
      dm[cc,] = c(this.noise, mean(clusters$Ar[I], na.rm=T), 
                  sum(clusters$Ar[I]>0.65, na.rm=T) / sum(!is.na(clusters$Ar[I])),
                  this.algorithm, this.size)
    }
  }
}
dm = dm[1:cc,]
dm$x = as.numeric(dm$x)
dm$y = as.numeric(dm$y)
dm$y2 = as.numeric(dm$y2)
tmp = sort(unique(clusters$size.factor))
dm$size.factor = tmp[as.numeric(dm$size.factor)]


ggplot(dm, aes(x=x, y=y2, color=size.factor)) + 
  geom_line(alpha=0.6, size=2) +  facet_wrap(~algorithm) +
  xlab("Interactome FPR") + ylab("Fraction of clusters\nreproducible (Ai>0.65)") + 
  coord_cartesian(ylim=c(0,1)) + theme_bw() 
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4A_01_v01.pdf"
ggsave(fn,width=10, height=3)

