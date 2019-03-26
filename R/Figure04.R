
require(ggplot2)
source("/Users/Mercy/Academics/Foster/ClusterExplore/R/functions.R")

load("/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr.Rda")
load("/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr_df_z.Rda")
if (0) {
  fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr.txt"
  clusters = as.data.frame(read_tsv(fn))
  clusters$noise_mag = as.numeric(clusters$noise_mag)
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
  
  # calculate Ai for each cluster in iter=1
  clusters$Ar = numeric(nrow(clusters))
  I.iter1 = which(clusters$iter==1)
  for (ii in 1:length(I.iter1)) {
    print(paste("progress =", ii, "/", length(I.iter1)))
    this.cluster = unlist(strsplit(clusters$cluster[I.iter1[ii]], ";"))
    this.algorithm = clusters$algorithm[I.iter1[ii]]
    this.fdr = clusters$noise_mag[I.iter1[ii]]
    
    # calculate Ai, complex reproducibility
    unqIters = 2:10
    Ar = numeric(length(unqIters))
    for (jj in 1:length(unqIters)) {
      I = clusters$algorithm %in% this.algorithm & clusters$noise_mag==this.fdr & clusters$iter==unqIters[jj]
      Ar[jj] = calcA(this.cluster, clusters$cluster[I])
    }
    
    clusters$Ar[I.iter1[ii]] = mean(Ar, na.rm=T)
  }
  save(clusters, file = "/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr.Rda")
  
  
  # calculate NULL Ai
  sizeRange = c(3,6,12)
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
  save(df.z, file = "/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr_df_z.Rda")
}



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
        clusters$size.factor%in%this.size & clusters$iter==1
      if (sum(I)==0) next
      dm[cc,] = c(this.noise, mean(clusters$Ar[I], na.rm=T), this.algorithm, this.size)
    }
  }
}
dm = dm[1:cc,]
dm$x = as.numeric(dm$x)
dm$y = as.numeric(dm$y)
tmp = sort(unique(clusters$size.factor))
dm$size.factor = tmp[as.numeric(dm$size.factor)]


Ar_null_n3 = mean(df.z$Ar[df.z$size==3],na.rm=T)
Ar_null_n6 = mean(df.z$Ar[df.z$size==6],na.rm=T)
Ar_null_n12 = mean(df.z$Ar[df.z$size==12],na.rm=T)


I.iter1 = clusters$iter==1
ggplot(clusters[I.iter1,], aes(x=noise_mag, y=Ar, color=size.factor)) + 
  geom_point(alpha=0.25) +  facet_wrap(~algorithm) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") + 
  geom_line(data=dm, aes(x=x, y=y,color=size.factor), size=2, alpha=.6) +
  coord_cartesian(ylim=c(0,1)) + theme_bw() +
  geom_hline(yintercept=Ar_null_n3, linetype="dashed", alpha=.65, colour="#F8766D") +
  geom_hline(yintercept=Ar_null_n6, linetype="dashed", alpha=.65, colour="#7CAE00")+
  geom_hline(yintercept=Ar_null_n12, linetype="dashed", alpha=.65, colour="#00BFC4")# + 
#annotate("text",x=.7,y=Ar_null_n3+.015, label="N=3", alpha=.75) +
#annotate("text",x=.7,y=Ar_null_n12+.015, label="N=12", alpha=.75)
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4A_01_v01.pdf"
ggsave(fn,width=10, height=3)


clusters$noise.factor = factor("fpr=0", levels = c("fpr=0", "fpr<=0.05","0.05<fpr<=0.50","0.50<fpr"))
clusters$noise.factor[clusters$noise_mag<=0.05] = as.factor("fpr<=0.05")
clusters$noise.factor[clusters$noise_mag==0] = as.factor("fpr=0")
clusters$noise.factor[clusters$noise_mag>0.05 & clusters$noise_mag<=0.5] = "0.05<fpr<=0.50"
clusters$noise.factor[clusters$noise_mag>=0.5] = "0.50<fpr"
ggplot(clusters[I.iter1,], aes(x=log10(size), y=Ar, color=noise.factor)) + 
  geom_jitter(alpha=0.15, width=0.015) + 
  geom_smooth() + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  coord_cartesian(ylim = c(-.001,1.001)) + theme_bw() + facet_wrap(~algorithm) +
  geom_hline(yintercept=Ar_null_n3, linetype="dashed", alpha=.65, colour="#F8766D") +
  geom_hline(yintercept=Ar_null_n6, linetype="dashed", alpha=.65, colour="#7CAE00")+
  geom_hline(yintercept=Ar_null_n12, linetype="dashed", alpha=.65, colour="#00BFC4")# + facet_grid(~algorithm)
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4B_01_v01.pdf"
ggsave(fn,width=10, height=4)


# zoom in on fdr=0.01
I = clusters$noise_mag %in% 0 | clusters$noise_mag %in% 0.01 | clusters$noise_mag %in% 0.02
ggplot(clusters[I&I.iter1,], aes(factor(noise_mag), y=Ar)) +  theme_bw() +
  geom_violin() + geom_jitter(width=.02,alpha=.05) + facet_grid(algorithm ~ size.factor) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") + 
  coord_cartesian(ylim=c(0,1))
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4X_01_v01.pdf"
ggsave(fn,width=10, height=6)



# hierarchical clusters
# reproduce Ji-size relationship for a toy dataset, then pick that dataset apart

# load chromatograms
fn = "/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/Input/Combined_replicates_2014_04_22_contaminates_removed_for_HvsL_scripts.csv"
chroms = as.data.frame(read_csv(fn))
tmp = names(chroms)
tmp[1] = "protid"
names(chroms) = tmp
# reduce chroms to replicate 3, good proteins
chroms = chroms[chroms$Replicate==2,]

# read proteasomal proteins
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/data/uniprot-proteasome+26s.tab"
uniprot = as.data.frame(read_tsv(fn))
I.26s = grepl("26S proteasome",uniprot$`Protein names`)
I.prot = (chroms$protid %in% uniprot$Entry[I.26s])

# read other complexes
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))
Icct = which(grepl("CCT",corum$ComplexName) & grepl("uman", corum$Organism))
prots.cct = unlist(strsplit(corum$`subunits(UniProt IDs)`[Icct], ";"))
Icct = chroms$protid %in% prots.cct
Isig = which(grepl("ignalos",corum$ComplexName) & grepl("uman", corum$Organism))
prots.sig = unlist(strsplit(corum$`subunits(UniProt IDs)`[Isig], ";"))
Isig = chroms$protid %in% prots.sig
Isplice = which(corum$ComplexName=="Spliceosome" & grepl("uman", corum$Organism))
prots.splice = unlist(strsplit(corum$`subunits(UniProt IDs)`[Isplice], ";"))
Isplice = chroms$protid %in% prots.splice

n.fill = rowSums(!is.na(chroms))
I.prot = unique(c(which(I.prot & n.fill>30), 
                  which(Icct & n.fill>30), 
                  which(Isig & n.fill>30), 
                  which(Isplice & n.fill>30),
                  #sample(which(n.fill==57), 1000)))
                  sample(which(n.fill>25), 50)))

# reduce chroms to proteasomal proteins and full proteins
chroms = chroms[I.prot,]

# fill NAs and remove v. high values
for (ii in 1:nrow(chroms)) {
  I = which(is.na(chroms[ii,]))
  tmp = chroms[ii,3:57]
  chroms[ii,I] = mean(as.numeric(tmp), na.rm=T)
  #chroms[ii,I] = 0
}

# hierarchical cluster
nclust = 25
chroms.mat = as.matrix(chroms[,3:57])
d0 = dist(chroms.mat)
hc0 = hclust(d0)
clusts0 = cutree(hc0,nclust)
I0 = order.hclust(hc0)
# make 10% noised copy
chroms.mat1 = chroms.mat * exp(matrix(0.1 * rnorm(nrow(chroms.mat) * ncol(chroms.mat)), 
                                       nrow = nrow(chroms.mat),
                                       ncol = ncol(chroms.mat)))
d1 = dist(chroms.mat1)
hc1 = hclust(d1)
clusts1 = cutree(hc1,nclust)
I1 = order.hclust(hc1)


# assemble complexes and calculate Ji, size
set0 = character(nclust)
set1 = character(nclust)
for (ii in 1:nclust) {
  I = which(clusts0==ii)
  set0[ii] = paste(I, collapse = ";")
  if (length(I)<3) set0[ii] = ""
  
  I = which(clusts1==ii)
  set1[ii] = paste(I, collapse = ";")
  if (length(I)<3) set1[ii] = ""
}
# remove clusters of size<3
set0 = set0[!set0==""]
set1 = set1[!set1==""]

df = data.frame(JJ = numeric(length(set1)), size = numeric(length(set1)), stringsAsFactors = F)
for (ii in 1:length(set1)) {
  df$JJ[ii] = calcA(set1[ii], set0)
  df$size[ii] = length(unlist(strsplit(set1[ii], ";")))
}

ggplot(df, aes(x=size, y= JJ)) + geom_point()

#heatmaply(scale(chroms.mat10), k_row = 15, Colv = F)

df = melt(chroms[I0,], id.vars=c("protid","Replicate"))
df$fraction = as.numeric(unlist(sapply(sapply(as.character(df$Var2), strsplit, " "),"[", 3)))
ggplot(df, aes(x=fraction, y=protid, fill=log(value))) + geom_tile()


