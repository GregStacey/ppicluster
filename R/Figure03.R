
source("functions.R")

# load chromatograms
fn = "../data/Combined_replicates_2014_04_22_contaminates_removed_for_HvsL_scripts.csv"
chroms = as.data.frame(read_csv(fn))
tmp = names(chroms)
tmp[1] = "protid"
names(chroms) = tmp

# load clusters
fn1 = "../data/clusters_wshuffle_moredata.txt"
fn2 = "../data/cluster3.txt"
data.c1 = as.data.frame(read_tsv(fn1))
data.c2 = as.data.frame(read_tsv(fn2))
data.c = rbind(data.c1, data.c2)
data.c$noise_mag = as.numeric(data.c$noise_mag)
# reduce to just experiment data
data.c = data.c[!data.c$data_type%in%"corum",]

# load interactomes
fn = "../data/interactomes_moredata.txt"
ints.c = as.data.frame(read_tsv(fn))
ints.c$ppi = paste(ints.c$protA, ints.c$protB, sep="-")

unqalgs = unique(data.c$algorithm)
unqdatasets = unique(data.c$data_type)
unqmags = sort(unique(data.c$noise_mag))
unqmags = unqmags[unqmags<=1]

# calculate Ai
fn = "../data/intJ_vs_clustJ_v03.Rda"
if (T){
  load(fn)
  df = df[!df$algorithm == "hierarchical",]
} else {
  nn = 10^6
  df = data.frame(dataset = character(nn), algorithm = character(nn), noise_mag=numeric(nn),
                  nclust0 = numeric(nn), nclust = numeric(nn), clust.size=numeric(nn), clustJ=numeric(nn),
                  nint0 = numeric(nn), nint = numeric(nn), intJ=numeric(nn),
                  stringsAsFactors = F)
  cc = 0
  for (ii in 1:length(unqdatasets)) {
    print(paste("   ", ii))
    i0 = ints.c$dataset%in%unqdatasets[ii] & ints.c$noise_mag==0
    this.nint0 = sum(i0)
    for (jj in 1:length(unqmags)) {
      print(unqmags[jj])
      i1 = ints.c$dataset%in%unqdatasets[ii] & ints.c$noise_mag==unqmags[jj]
      this.intJ = length(intersect(ints.c$ppi[i0], ints.c$ppi[i1])) / 
        length(unique(c(ints.c$ppi[i0],ints.c$ppi[i1])))
      this.nint = sum(i1)
      for (kk in 1:length(unqalgs)) {
        I0 = data.c$data_type%in%unqdatasets[ii] & data.c$noise_mag==0 & data.c$algorithm%in%unqalgs[kk]
        I1 = data.c$data_type%in%unqdatasets[ii] & data.c$noise_mag==unqmags[jj] & data.c$algorithm%in%unqalgs[kk]
        if (sum(I1)==0) next
        
        cluster0 = data.c$cluster[I0]
        cluster = data.c$cluster[I1]
        J.i = numeric(length(cluster))
        for (mm in 1:length(cluster)) {
          J.i[mm] = calcA(cluster[mm], cluster0)
          
          cc = cc+1
          df$dataset[cc] = unqdatasets[ii]
          df$algorithm[cc] = unqalgs[kk]
          df$noise_mag[cc] = unqmags[jj]
          df$nclust0[cc] = length(cluster0)
          df$nclust[cc] = length(cluster)
          df$clustJ[cc] = J.i[mm]
          df$nint0[cc] = this.nint0
          df$nint[cc] = this.nint
          df$intJ[cc] = this.intJ
          df$clust.size[cc] = length(unlist(strsplit(cluster[mm], ";")))
        }
        
      }
    }
  }
  df = df[1:cc,]
  
  save("df", file = fn)
}

df$experiment = as.numeric(unlist(sapply(sapply(df$dataset, strsplit, "_"), "[", 2)))
groups = list(c(1,2,3,4,5,6,7,8),c(9,10,11,12,13,14),c(15,16,17,18,19,20),c(21,22,23,24,25,26,27,28))
df$experiment[df$experiment%in%groups[[1]]] = "DS1"
df$experiment[df$experiment%in%groups[[2]]] = "DS2"
df$experiment[df$experiment%in%groups[[3]]] = "DS3"
df$experiment[df$experiment%in%groups[[4]]] = "DS4"
df$algorithm[df$algorithm=="co"] = "CO"
df$algorithm[df$algorithm=="co_mcl"] = "CO+MCL"
df$algorithm[df$algorithm=="mcl"] = "MCL"
df$algorithm[df$algorithm=="pam"] = "k-Med"


# get averages for figures
nn = 10^3
dm = data.frame(x = numeric(nn), 
                y = numeric(nn),
                intJ = numeric(nn),
                nclust = numeric(nn),
                algorithm = character(nn),
                experiment = character(nn),
                stringsAsFactors = F)
cc = 0
for (ii in 1:length(unique(df$noise_mag))) {
  this.noise = sort(unique(df$noise_mag))[ii]
  for (jj in 1:length(unique(df$algorithm))) {
    this.algorithm = sort(unique(df$algorithm))[jj]
    for (kk in 1:length(unique(df$experiment))) {
      this.experiment = sort(unique(df$experiment))[kk]
      cc = cc+1
      I = df$noise_mag==this.noise & 
        df$algorithm%in%this.algorithm & 
        df$experiment%in%this.experiment
      dm[cc,] = c(this.noise, mean(df$clustJ[I], na.rm=T), 
                  mean(df$intJ[I], na.rm=T), sum(I), this.algorithm, this.experiment)
    }
  }
}
dm = dm[1:cc,]
dm$x = as.numeric(dm$x)
dm$y = as.numeric(dm$y)
dm$intJ = as.numeric(dm$intJ)
dm$nclust = as.numeric(dm$nclust)




# A. Proteasome chromatograms
# read proteasomal proteins
fn = "../data/uniprot-proteasome+26s.tab"
uniprot = as.data.frame(read_tsv(fn))
I.26s = grepl("26S proteasome",uniprot$`Protein names`)
I = chroms$protid %in% uniprot$Entry[I.26s]
df.prot = melt(chroms[I,], id.vars=c("protid","Replicate"))
df.prot$variable = gsub(" ", "..", df.prot$variable)
df.prot$fraction = as.numeric(unlist(sapply(sapply(df.prot$variable, strsplit, "..", fixed=T), "[", 3)))
df.prot = df.prot[,-3]
df.prot$value2 = df.prot$value * exp(rnorm(nrow(df.prot)) * 0.25)
df.prot$value3 = df.prot$value * exp(rnorm(nrow(df.prot)) * 0.5)
df.prot = melt(df.prot, id.vars=c("protid", "Replicate", "fraction"))
df.prot$variable = as.character(df.prot$variable)
df.prot$variable[df.prot$variable %in% "value"] = "No noise"
df.prot$variable[df.prot$variable %in% "value2"] = "25% noise"
df.prot$variable[df.prot$variable %in% "value3"] = "50% noise"
df.prot$variable = factor(df.prot$variable, levels = c("No noise", "25% noise", "50% noise"))
df.prot = df.prot[df.prot$Replicate==1,]

# calculate average R^2
unqprots = unique(df.prot$protid)
unqnoise = unique(df.prot$variable)
for (ii in 1:length(unqnoise)) {
  Rmat = matrix(nrow = length(unqprots), ncol = length(unqprots))
  for (jj in 1:length(unqprots)) {
    ia = which(df.prot$protid==unqprots[jj] & df.prot$variable==unqnoise[ii])
    for (kk in 1:length(unqprots)) {
      if (jj>=kk) next
      ib = which(df.prot$protid==unqprots[kk] & df.prot$variable==unqnoise[ii])
      I = !is.na(df.prot$value[ia]) & !is.na(df.prot$value[ib])
      if (sum(I)<3) next
      Rmat[jj,kk] = cor.test(log(df.prot$value[ia]), log(df.prot$value[ib]), na.rm=T)$estimate
    }
  }
  Rmat[Rmat==0] = NA
  print(paste(unqnoise[ii], mean(Rmat, na.rm=T)^2))
}


ggplot(df.prot[df.prot$Replicate==1,], aes(x=fraction,y=log(value),group=protid)) + geom_line(alpha=.35) +
  facet_wrap(~variable) + theme_bw() + xlab("Fraction") + ylab("log(Protein amount)")
fn = "../figures/fig_3A_v01.pdf"
ggsave(fn,width=10, height=3)




# B. Ai vs chromnoise
ggplot(df, aes(x=noise_mag*100, y=clustJ, color=experiment)) + 
  geom_point(alpha=0.01) +  facet_grid(~algorithm) +
  ylab("Similarity to un-noised (Ji)") + xlab("Noise magnitude, %") + 
  geom_line(data=dm, aes(x=x*100, y=y,color=experiment), size=2, alpha=.6) +
  theme_bw() + xlim(0,50) + theme(legend.position = "none") +
  scale_color_brewer(palette="Set1")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_3B_v03.pdf"
ggsave(fn,width=10, height=3)
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_3B_v03.png"
ggsave(fn,width=10, height=3)



# C. compare interactome changes to cluster changes
I = which(df$nint>1000 & df$nint0>1000 & df$nclust>10)
I2 = dm$experiment%in%"group1" & dm$intJ<0.1 & dm$y>0.5
dm$y[I2] = NA
ggplot(df[I,], aes(x=intJ, y=clustJ, color=experiment)) + 
  geom_point(alpha=0.05) +  facet_grid(~algorithm) +
  ylab("Cluster similarity, J") + xlab("Interactome similarity, Jaccard") + 
  geom_line(data=dm[dm$nclust>10,], aes(x=intJ, y=y,color=experiment), size=2, alpha=.6) +
  theme_bw()
fn = "../figures/fig_3C_v02.pdf"
ggsave(fn,width=10, height=3)
