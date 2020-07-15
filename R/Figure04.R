
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

# add walktrap
fn.walk = "../data/data_c_walktrap.txt"
if (T) {
  load(fn.walk)
} else {
  # walktrap cluster
  # add 10 iterations of PAM and walktrap
  data.c.add = data.frame(data_type = character(10^5), noise_type = character(10^5),
                          noise_mag = numeric(10^5), algorithm = character(10^5),
                          cluster = character(10^5), stringsAsFactors = F)
  cc = 0
  for (ii in 1:length(unqmags)) {
    for (jj in 1:length(unqdatasets)) {
      print(paste("walktrap clustering dataset",unqdatasets[jj],"at noise=", unqmags[ii]))
      
      # get shuffled network
      I = ints.c$noise_mag==unqmags[ii] & ints.c$dataset==unqdatasets[jj] & !ints.c$protA==ints.c$protB
      if (sum(I)<190) next
      these.ints = ints.c[I,c("protA", "protB")]
      ints.shuffle = shufflecorum(these.ints, unqmags[ii])
      
      # walktrap
      graph.object = graph_from_edgelist(as.matrix(ints.shuffle), directed = F)
      walk.cluster = walktrap.community(graph.object)
      for (kk in 1:length(walk.cluster)) {
        if (length(walk.cluster[[kk]]) < 3) next
        cc = cc+1
        data.c.add$data_type[cc] = unqdatasets[jj]
        data.c.add$noise_type[cc] = "chrom"
        data.c.add$noise_mag[cc] = unqmags[ii]
        data.c.add$algorithm[cc] = "walk"
        data.c.add$cluster[cc] = paste(walk.cluster[[kk]], collapse=";")
      }
    }
  }
  data.c.add = data.c.add[1:cc,]
}
data.c = rbind(data.c, data.c.add)


# calculate Ji
fn = "../data/intJ_vs_clustJ_v03.Rda"
fn.walk = "../data/intJ_vs_clustJ_v03walk.Rda"
if (T){
  load(fn.walk)
  tmp = df
  load(fn)
  df = rbind(df, tmp)
  #df = df[!df$algorithm == "hierarchical",]
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
df = df[, -11]

# add louvain, and leiden
tmp = as.data.frame(read_tsv("../data/data.c.add_cofracmcp_summary.txt"))
#tmp = tmp[,6:15]
#names(tmp) = names(df)
df = rbind(df, tmp)
# add mcode
tmp = as.data.frame(read_tsv("../data/data.c.add_cofracmcp_summary_mcode.txt"))
df = rbind(df, tmp)


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
df$algorithm[df$algorithm=="mcode"] = "MCODE"
df$algorithm[df$algorithm=="walk"] = "walktrap"
df$algorithm[df$algorithm=="hierarchical"] = "Hierarchical"
df$algorithm[df$algorithm=="louvain"] = "Louvain"
df$algorithm[df$algorithm=="leiden"] = "Leiden"

#write_tsv(df, path="../data/interactomes_moredata.txt")


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


# calculate R^2 as function of chrom_noise
noise.range = seq(from=0, to=1.25, by=0.01)
chrom.mat0 = as.matrix(chroms[,3:57])
good.data = !is.na(chrom.mat0)
n.data = rowSums(good.data)
n.data[sample(length(n.data), round(length(n.data)*.9))] = 0
nn = sum(n.data>5)
chrom.mat0 = chrom.mat0[n.data>5,]
good.data = good.data[n.data>5,]
df.chromcor = data.frame(chromnoise = numeric(10^6), 
                         rr = numeric(10^6),
                         avg.rr = numeric(10^6))
for (ii in 1:length(noise.range)) {
  print(noise.range[ii])
  noise.mat = matrix(rnorm(nrow(chrom.mat0) * ncol(chrom.mat0)), 
                     nrow = nrow(chrom.mat0), ncol = ncol(chrom.mat0))
  chrom.mat = chrom.mat0 * exp(noise.mat * noise.range[ii])
  
  RR = numeric(nrow(chrom.mat))
  for (jj in 1:nrow(chrom.mat)) {
    RR[jj] = cor(log(chrom.mat0[jj,good.data[jj,]]), 
                 log(chrom.mat[jj,good.data[jj,]]))
  }
  RR = RR ^ 2
  
  I = ((ii-1)*nn+1) : (ii*nn)
  df.chromcor$rr[I] = RR
  df.chromcor$avg.rr[I] = mean(RR, na.rm=T)
  df.chromcor$chromnoise[I] = noise.range[ii]
}
df.chromcor = df.chromcor[1:max(I),]




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
  RR = rep(NA, length(unqprots))
  for (jj in 1:length(unqprots)) {
    I0 = which(df.prot$protid == unqprots[jj] & df.prot$variable==unqnoise[1])
    I1 = which(df.prot$protid == unqprots[jj] & df.prot$variable==unqnoise[ii])
    
    ia = !is.na(df.prot$value[I0])
    if (sum(ia)<3) next
    RR[jj] = cor.test(log(df.prot$value[I0[ia]]), log(df.prot$value[I1[ia]]))$estimate
  }
  print(paste(unqnoise[ii], mean(RR^2, na.rm=T)))
}


ggplot(df.prot[df.prot$Replicate==1,], aes(x=fraction,y=log(value),group=protid)) + geom_line(alpha=.35) +
  facet_wrap(~variable) + theme_bw() + xlab("Fraction") + ylab("log(Protein amount)")
fn = "../figures/fig_3A_v01.pdf"
ggsave(fn,width=10, height=3)



# B. Ji vs chromnoise
ggplot(df, aes(x=noise_mag*100, y=clustJ, color=experiment)) + 
  geom_point(alpha=0.01) +  facet_wrap(~algorithm, nrow = 2) +
  ylab("Similarity to un-noised (Ji)") + xlab("Co-fractionation noise, %") + 
  geom_line(data=dm, aes(x=x*100, y=y,color=experiment), size=2, alpha=.6) +
  theme_bw() + xlim(0,50) + theme(legend.position = "none") + scale_colour_grey()
#scale_color_brewer(palette="Set1")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_3B_v05.pdf"
ggsave(fn,width=10, height=6)
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_3B_v05.png"
ggsave(fn,width=10, height=6)



# C. chromatogram
ggplot(df.chromcor, aes(x=chromnoise*100, y=rr)) + 
  geom_point(alpha=.002) + geom_line(aes(x=chromnoise*100, y=avg.rr)) +
  xlab("Noise magnitude, %") + ylab("Noised vs un-noised chromatogram\ncorrelation (Pearson R)") +
  theme_bw()
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_3C_v03.pdf"
ggsave(fn,width=3.4, height=3)
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_3C_v03.png"
ggsave(fn,width=3.4, height=3)



# D. co-interactome probability vs chromnoise

# make coint prob
unqsets = unique(ints.c$dataset)
unqmags = unique(ints.c$noise_mag)
nn = length(unqsets) * length(unqmags)
iterMax = 100
df.intJ = data.frame(dataset = character(nn),
                     noise_mag = numeric(nn),
                     avg.coint = numeric(nn), 
                     avg.intJ = numeric(nn), 
                     coint = numeric(nn), 
                     intJ = numeric(nn), stringsAsFactors = F)
df.intJ[df.intJ==0] = NA
cc = 0
for (ii in 1:length(unqsets)) {
  I0 = ints.c$dataset==unqsets[ii] & ints.c$noise_mag==0
  ints0 = ints.c$ppi[I0]
  if (length(ints0)<250) next
  for (jj in 1:length(unqmags)) {
    I = ints.c$dataset==unqsets[ii] & ints.c$noise_mag==unqmags[jj]
    ints = ints.c$ppi[I]
    if (length(ints)<250) next
    cc = cc+1
    df.intJ$coint[cc] = sum(ints0 %in% ints) / length(ints0)
    df.intJ$intJ[cc] = length(intersect(ints0, ints)) / length(unique(c(ints0, ints)))
    df.intJ$dataset[cc] = unqsets[ii]
    df.intJ$noise_mag[cc] = unqmags[jj]
  }
}
for (ii in 1:length(unqmags)) {
  I = df.intJ$noise_mag==unqmags[ii]
  df.intJ$avg.coint[I] = mean(df.intJ$coint[I], na.rm=T)
  df.intJ$avg.intJ[I] = mean(df.intJ$intJ[I], na.rm=T)
}

ggplot(df.intJ, aes(x=noise_mag*100, y=coint)) + geom_point(alpha=.1) +
  geom_line(aes(x=noise_mag*100, y=avg.coint)) +
  xlab("Noise magnitude, %") + ylab("Interactome Jaccard") +
  theme_bw()
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_3D_v05.pdf"
ggsave(fn,width=3.4, height=3)
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_3D_v05.png"
ggsave(fn,width=3.4, height=3)

