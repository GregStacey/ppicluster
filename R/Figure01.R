##### ------------------------------------------- #####
# 0. initialize

set.seed(4)

# make your own melt function
meltmat = function(mat, id.vars) {
  nmes = names(mat)
  I = !nmes %in% id.vars
  mat.to.melt = mat[,I]
  
  nn = nrow(mat.to.melt)
  nn2 = nn^2
  
  df = data.frame(x=numeric(nn2), y=numeric(nn2), value=numeric(nn2), proteins=numeric(nn2))
  cc = 0
  for (ii in 1:nrow(mat.to.melt)) {
    for (jj in 1:ncol(mat.to.melt)) {
      cc = cc+1
      df$x[cc] = ii
      df$y[cc] = jj
      df$value[cc] = mat.to.melt[ii,jj]
      df$proteins[cc] = mat$proteins[ii]
    }
  }
  
  return(df)
}

# get corum
fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))
corum = corum[corum$Organism%in%"Human",]
corum$size = unlist(lapply(lapply(lapply(corum$`subunits(UniProt IDs)`, strsplit, ";"), unlist), length))


##### ------------------------------------------- #####
# plot example ppi network
I = sample(which(corum$size>3), 6)
prots = unlist(lapply(corum$`subunits(UniProt IDs)`[I], strsplit, ";"))
conmat = matrix(nrow=length(prots), ncol=length(prots))
for (ii in 1:length(I)) {
  this.prots = unlist(strsplit(corum$`subunits(UniProt IDs)`[I[ii]], ";"))
  this.I = which(prots %in% this.prots)
  for (jj in 1:length(this.I)) {
    for (kk in 1:length(this.I)) {
      if (jj==kk) next
      conmat[this.I[jj], this.I[kk]] = 1
    }
  }
}
conmat.df = as.data.frame(conmat)
names(conmat.df) = as.character(1:length(prots))
conmat.df$proteins = as.character(1:length(prots))
df = meltmat(conmat.df, id.vars = c("proteins"))
df$value[is.na(df$value)] = 0

ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme_bw() + blank_theme + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_1A_1_v01.png",
       width=4, height=4)




##### ------------------------------------------- #####
# plot noised ppi network
precRange = c(0.9, 0.5)
I1 = which(conmat==1)
I0 = which(conmat==0)
mm = nrow(noisemat)
for (jj in 1:length(precRange)) {
  # just noise
  prec = precRange[jj]
  I.out = I1[sample(1:length(I1), round(length(I1) * (1-prec) / 2))]
  noisemat = matrix(0, nrow(conmat), ncol(conmat))
  for (ii in 1:length(I.out)) {
    I.in = I0[sample(1:length(I0), 1)]
    I.in.a = ((I.in-1) %% mm) + 1
    I.in.b = floor((I.in-1) / mm) + 1
    I.out.a = ((I.out[ii]-1) %% mm) + 1
    I.out.b = floor((I.out[ii]-1) / mm) + 1
    
    # upper triangular
    noisemat[I.out.a, I.out.b] = -1
    noisemat[I.in.a, I.in.b] = 1

    # lower triangular
    noisemat[I.out.b, I.out.a] = -1
    noisemat[I.in.b, I.in.a] = 1
  }
  noisemat.df = as.data.frame(noisemat)
  names(noisemat.df) = as.character(1:length(prots))
  noisemat.df$proteins = as.character(1:length(prots))
  df.noise = meltmat(noisemat.df, id.vars = c("proteins"))
  df.noise$value[is.na(df.noise$value)] = 0
  
  ggplot(df.noise, aes(x, y)) +
    geom_tile(aes(fill = abs(value)), color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") + 
    theme_bw() + blank_theme + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  fn = paste("/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_1A_", (jj-1)*2+2, "_v01.png", sep="")
  ggsave(fn,width=4, height=4)
  
  
  # network+noise
  noisemat[is.na(noisemat)] = 0
  conmat[is.na(conmat)] = 0
  fullmat = noisemat + conmat
  fullmat[fullmat==2] = 0
  fullmat.df = as.data.frame(fullmat)
  names(fullmat.df) = as.character(1:length(prots))
  fullmat.df$proteins = as.character(1:length(prots))
  df.full = meltmat(fullmat.df, id.vars = c("proteins"))
  df.full$value[is.na(df.full$value)] = 0
  
  ggplot(df.full, aes(x, y)) +
    geom_tile(aes(fill = value), color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") + 
    theme_bw() + blank_theme + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  fn = paste("/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_1A_", (jj-1)*2+3, "_v01.png", sep="")
  ggsave(fn,width=4, height=4)
}




##### ------------------------------------------- #####
# Simple Ai vs FDR

# get clusters
fn.save = "/Users/Mercy/Academics/Foster/ClusterExplore/data/figure01_panel01.Rda"
if (TRUE){
  load(fn.save)
} else {
  # define functions
  source("/Users/Mercy/Academics/Foster/ClusterExplore/R/functions.R")
  
  fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr.txt"
  data = read_tsv(fn)
  clusters = data[data$iter==1 & data$algorithm%in%"mcl",]
  
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
    for (jj in 1:length(unqIters)) {
      I = data$algorithm %in% this.algorithm & data$noise_mag==this.fdr & data$iter==unqIters[jj]
      Ar[jj] = calcA(this.cluster, data$cluster[I])
    }
    clusters$Ar[I.order[ii]] = mean(Ar, na.rm=T)
  }
  
  # get averages for figures
  nn = 10^3
  dm = data.frame(x = numeric(nn), 
                  y = numeric(nn),
                  algorithm = character(nn),stringsAsFactors = F)
  cc = 0
  for (ii in 1:length(unique(clusters$noise_mag))) {
    this.noise = sort(unique(clusters$noise_mag))[ii]
    for (jj in 1:length(unique(clusters$algorithm))) {
      this.algorithm = sort(unique(clusters$algorithm))[jj]
      
      cc = cc+1
      I = clusters$noise_mag==this.noise & 
        clusters$algorithm%in%this.algorithm 
      dm[cc,] = c(this.noise, mean(clusters$Ar[I], na.rm=T), this.algorithm, this.size)
      
    }
  }
  dm = dm[1:cc,]
  dm$x = as.numeric(dm$x)
  dm$y = as.numeric(dm$y)
  
  save(clusters, dm, file=fn.save)
}

ggplot(clusters, aes(x=noise_mag, y=Ar)) + geom_point(alpha=0.04,color='steelblue') +
  ylab("Cluster-wise similarity, Ji") + xlab("Network FPR") + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_line(data=dm, aes(x=x, y=y), size=1, alpha=.9,color='steelblue') + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()) 
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_1A_6_v01.png"
ggsave(fn,width=5, height=2.5)


