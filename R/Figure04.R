
require(ggplot2)
source("functions.R")

#fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/data/fig2B_v03.Rda"
#load(fn) # sim, Ji
fn = "../data/clusters_full_netw2.txt"
Ji = as.data.frame(read_tsv(fn))

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


Ar_null_n3 = mean(df.z$Ar[df.z$size==3],na.rm=T)
Ar_null_n6 = mean(df.z$Ar[df.z$size==6],na.rm=T)
Ar_null_n12 = mean(df.z$Ar[df.z$size==12],na.rm=T)


ggplot(Ji[Ji$iter>1,], aes(x=noise_mag, y=Ji2, color=size.factor)) + 
  geom_point(alpha=0.05) + facet_wrap(~algorithm) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") + 
  geom_line(data=dm, size=2, alpha=.6) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,0.5)) + theme_bw() +
  geom_hline(yintercept=Ar_null_n3, linetype="dashed", alpha=.65, colour="#F8766D") +
  geom_hline(yintercept=Ar_null_n6, linetype="dashed", alpha=.65, colour="#7CAE00")+
  geom_hline(yintercept=Ar_null_n12, linetype="dashed", alpha=.65, colour="#00BFC4")
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4A_01_v01.pdf"
#ggsave(fn,width=10, height=3)


ggplot(Ji[Ji$cluster.size<150 & Ji$iter>1,], aes(x=log10(cluster.size), y=Ji2)) + 
  geom_jitter(alpha=0.05, width=0.03, height=0.03) + facet_grid(noise.factor~algorithm) +
  geom_smooth() + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  coord_cartesian(ylim = c(-.001,1.001)) + theme_bw() + 
  geom_hline(yintercept=Ar_null_n3, linetype="dashed", alpha=.65, colour="#F8766D") +
  geom_hline(yintercept=Ar_null_n6, linetype="dashed", alpha=.65, colour="#7CAE00")+
  geom_hline(yintercept=Ar_null_n12, linetype="dashed", alpha=.65, colour="#00BFC4")
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4B_01_v01.pdf"
#ggsave(fn,width=10, height=4)


# zoom in on fdr=0.01
I = Ji$noise_mag %in% 0 | Ji$noise_mag %in% 0.01 | Ji$noise_mag %in% 0.02
ggplot(Ji[I,], aes(factor(noise_mag), y=Ji1)) +  theme_bw() +
  geom_violin() + geom_jitter(width=.02,alpha=.05) + facet_grid(algorithm ~ size.factor) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") + 
  coord_cartesian(ylim=c(0,1))
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4X_01_v01.pdf"
#ggsave(fn,width=10, height=6)


# bar graph of smallest 20% of clusters vs top 80% (or something)
tmp = Ji
tmp$size.split = factor(nrow(tmp), ordered = T)
levels(tmp$size.split) = c("low", "high")
tmp$size.split[tmp$cluster.size<=15] = "low"
tmp$size.split[tmp$cluster.size>15] = "high"
#tmp$size.split[tmp$size.prctile<=0.05] = "low"
#tmp$size.split[tmp$size.prctile>0.05] = "high"
ggplot(tmp, aes(y=Ji2, fill=size.split, x=algorithm)) + geom_boxplot()


# size distributions
ggplot(Ji[Ji$iter>1,], aes(y=log10(cluster.size), x=algorithm)) + geom_violin()


# stats
I = Ji$iter>1
glm(Ji2 ~ noise_mag + algorithm + cluster.size, gaussian, Ji[I,])








ggplot(Ji[Ji$algorithm=="hierarchical" & Ji$cluster.size<500,], 
       aes(x=log10(jitter(cluster.size,factor=0.1)), y=Ji1)) + 
  geom_point(alpha=0.05) + facet_grid(noise.factor~iter) + 
  geom_smooth(method=lm, formula = y~x) +
  coord_cartesian(ylim=c(0,1))






# hierarchical clusters
# reproduce Ji-size relationship for a toy dataset, then pick that dataset apart

# load chromatograms
fn = "/Users/Mercy/Academics/Foster/NickCodeData/
GregPCP-SILAC/Input/Combined_replicates_2014_04_22_contaminates_removed_for_HvsL_scripts.csv"
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


