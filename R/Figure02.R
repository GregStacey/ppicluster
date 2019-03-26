
require(ggplot2)
require(viridis)
require(scatterpie)
source("/Users/Mercy/Academics/Foster/ClusterExplore/R/functions.R")

blank_theme = theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(), legend.position="none")

fnsave = "/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_wshuffle_coexp.Rda"
load(fnsave)
clusters$noise_mag = as.numeric(clusters$noise_mag)

fnsave = "/Users/Mercy/Academics/Foster/ClusterExplore/data/cluster3.txt"
clusters2 = as.data.frame(read_tsv(fnsave))
clusters2$noise_mag = as.numeric(clusters2$noise_mag)
clusters2 = clusters2[clusters2$data_type=="corum",]

clusters = rbind(clusters[,1:5], clusters2)



# 2A ##### ------------------------------------------- #####
# example of noise-adding where
#   i) clusters are obviously changing
#   ii)  measures are different

I = clusters$data_type%in%"corum" & clusters$algorithm%in%"mcl" & clusters$noise_type%in%"network_shuffle"
data = clusters[I,]
set0 = data$cluster[data$noise_mag==0]

unqnoise = unique(data$noise_mag)
for (uu in 1:length(unqnoise)) {
  set1 = data$cluster[data$noise_mag==unqnoise[uu]]
  
  pci = matrix(nrow=length(set1), ncol=length(set0))
  pci.sort = matrix(nrow=length(set1), ncol=length(set0))
  for (ii in 1:length(set1)) {
    
    # find best-matching cluster
    cluster1 = unlist(strsplit(set1[ii], ";"))
    JJ = numeric(length(set0))
    for (jj in 1:length(set0)) {
      cluster0 = unlist(strsplit(set0[jj], ";"))
      nn.jacc = length(unique(c(cluster0, cluster1)))
      JJ[jj] = length(intersect(cluster0, cluster1)) / nn.jacc 
    }
    I.max = which.max(JJ)[1]
    cluster0 = unlist(strsplit(set0[I.max], ";"))
    nn.jacc = length(unique(c(cluster0, cluster1)))
    nn.intersect = length(intersect(cluster0, cluster1))
    
    # where are the rest of the proteins?
    non.overlapping = cluster1[!cluster1 %in% cluster0]
    N.overlap = numeric(length(set0))
    if (length(non.overlapping>0)) {
      for (jj in 1:length(non.overlapping)) {
        I = which(grepl(non.overlapping[jj], set0))
        if (length(I)>=1) {
          I = sample(I,1)
          N.overlap[I] = N.overlap[I] + 1
        }
      }
    }
    
    # store in pci
    pci[ii,] = N.overlap
    pci[ii,I.max] = JJ[I.max] * nn.jacc
    
    # order
    pci.sort[ii,] = sort(pci[ii,], decreasing = T)
    # inject columns: 
    pci.sort[ii,(-1:0)+ncol(pci.sort)] = pci.sort[ii,1] + 1
    pci.sort[ii,] = sort(pci.sort[ii,], decreasing = T)
    #   column j=1: unmatched fraction of cluster 0
    pci.sort[ii,1] = length(cluster0) - nn.intersect
    #   column j=2: unmatched fraction of cluster 1
    pci.sort[ii,2] = length(cluster1) - sum(pci.sort[ii,3:ncol(pci.sort)])
    
    # normalize
    pci[ii,] = pci[ii,] / nn.jacc
    pci.sort[ii,] = pci.sort[ii,] / nn.jacc
  }
  pci = as.data.frame(pci)
  pci$cluster = 1:nrow(pci)
  pci$size = unlist(lapply(lapply(lapply(set1, strsplit, ";"), "[[", 1), length))
  pci.sort = as.data.frame(pci.sort)
  pci$cluster = 1:nrow(pci.sort)
  pci.sort$size = unlist(lapply(lapply(lapply(set1, strsplit, ";"), "[[", 1), length))

  # tile x and y
  # assume that cluster with sqrt(size) needs a box of width=m*sqrt(size)
  gridx = 250
  df = pci.sort[,]
  df = df[order(df$size, decreasing = T),]
  mm = 3
  df$tile.size = round(mm * sqrt(df$size)/5)*5 * .85 + 2
  df$y = numeric(nrow(df))
  df$x = cumsum(df$tile.size) - df$tile.size[1]*.9
  df$x = (df$x %% gridx)
  y0 = 0
  I0 = 1
  for (ii in 2:nrow(df)) {
    if (abs(df$x[ii-1] - df$x[ii]) > gridx*.65) { 
      y0 = y0 + max(df$tile.size[I0:ii])/2 + 5
      I0 = ii
    }
    df$y[ii] = y0
  }
  df = df[,c(1:25,ncol(df) + seq(from=-3, to=0, by=1))]
  df[df==0] = 10^-4
  
  # plot
  ggplot() + geom_scatterpie(aes(x=x, y=y, r=sqrt(size)), data=df, cols=names(df)[1:25]) + 
    coord_fixed() + theme_void() + theme(legend.position="none")
  sf = paste("/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_", uu, ".pdf", sep="")
  ggsave(sf, width=10, height=5)
  sf = paste("/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_", uu, ".png", sep="")
  ggsave(sf, width=10, height=5, dpi=1000)
  
  # zoom in
  ggplot() + geom_scatterpie(aes(x=x, y=y, r=sqrt(size)), data=df, cols=names(df)[1:25]) + 
    coord_fixed() + theme_void() + theme(legend.position="none") + coord_cartesian(xlim = c(80, 130))
  sf = paste("/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_zoom_", uu, ".pdf", sep="")
  ggsave(sf, width=2, height=5)
  sf = paste("/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_zoom_", uu, ".png", sep="")
  ggsave(sf, width=2, height=5, dpi=1000)
}




# 2B ##### ------------------------------------------- #####
# line plots of A_i vs noise

sf = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/data/fig2B_v03.Rda"
if (F){
  load(sf)
} else {
  I = clusters$data_type%in%"corum" & clusters$noise_type%in%"network_shuffle"
  data = clusters[I,]
  unqalg = unique(data$algorithm)
  unqnoise = unique(data$noise_mag)
  sim = data.frame(measure = character(10^4),
                   algorithm = character(10^4),
                   noise_mag = character(10^4),
                   value = numeric(10^4), stringsAsFactors = F)
  Ji = data.frame(ji = numeric(10^6), algorithm = character(10^6), 
                  cluster.size = numeric(10^6), noise_mag=numeric(10^6),
                  stringsAsFactors = F)
  cc = 0
  cc2 = 0
  for (ii in 1:length(unqnoise)) {
    print(unqnoise[ii])
    for (jj in 3:length(unqalg)) {
      print(paste("  ",unqalg[jj]))
      I0 = data$noise_mag == 0 & data$algorithm%in%unqalg[jj]
      I1 = data$noise_mag == unqnoise[ii] & data$algorithm%in%unqalg[jj]
      set0 = data$cluster[I0]
      set1 = data$cluster[I1]
      if (sum(I1)==0 | sum(I0)==0) next
      
      # geomacc
      cc = cc+1
      sim$measure[cc] = "ga"
      sim$algorithm[cc] = unqalg[jj]
      sim$noise_mag[cc] = unqnoise[ii]
      sim$value[cc] = geomacc(set0, set1)
      
      # matching ratio
      cc = cc+1
      sim$measure[cc] = "mmr"
      sim$algorithm[cc] = unqalg[jj]
      sim$noise_mag[cc] = unqnoise[ii]
      sim$value[cc] = matchingratio(set0, set1)
      
      # nmi
      
      # Ai
      cc = cc+1
      sim$measure[cc] = "ai"
      sim$algorithm[cc] = unqalg[jj]
      sim$noise_mag[cc] = unqnoise[ii]
      tmp = numeric(length(set1))
      for (kk in 1:length(tmp)) {
        cc2 = cc2+1
        tmp[kk] = calcA(set1[kk], set0)
        Ji$ji[cc2] = tmp[kk]
        Ji$algorithm[cc2] = unqalg[jj]
        Ji$noise_mag[cc2] = unqnoise[ii]
        Ji$cluster.size[cc2] = length(unlist(strsplit(set1[kk], ";")))
      }
      sim$value[cc] = mean(tmp)
    }
  }
  sim = sim[1:cc,]
  sim$noise_mag = as.numeric(sim$noise_mag)
  sim$value[sim$value==0] = NA
  Ji = Ji[1:cc2,]
  Ji$noise_mag = as.numeric(Ji$noise_mag)

  save(list=c("sim","Ji"), file = sf)
}

# TO DO:
# add Ai scatter, e.g. using `clusters` from Figure03
# clusters2 = clusters
# names(clusters2) = c("iter", "noise_mag", "algorithm", "cluster" ,"value", "size", "size.factor", "noise.factor")
# clusters2$measure = "ai"
#   + geom_point(data=clusters2, alpha = 0.05)

Ji$measure = Ji$algorithm
Ji$value = Ji$ji
ggplot(sim, aes(x=noise_mag, y=value, color=measure)) + geom_line() + 
  facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
  geom_point(data = Ji, alpha=0.015, color="red")
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2B_v03.pdf"
ggsave(fn,width=10, height=3)



# 2C ##### ------------------------------------------- #####
# violin plots of A_i vs noise


sf = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/data/fig2C_v03.Rda"
if (T){
  load(sf)
} else {
  I = clusters$data_type%in%"corum" & clusters$noise_type%in%"network_shuffle" & clusters$noise_mag<=0.02
  data = clusters[I,]
  unqalg = unique(data$algorithm)
  unqnoise = c(0, 0.01, 0.02)
  ai = data.frame(algorithm = character(10^4),
                  noise_mag = character(10^4),
                  ai = numeric(10^4), stringsAsFactors = F)
  ai$ai[ai$ai==0] = NA
  cc = 0
  for (ii in 1:length(unqnoise)) {
    print(unqnoise[ii])
    for (jj in 1:length(unqalg)) {
      I0 = data$noise_mag == 0 & data$algorithm%in%unqalg[jj]
      I1 = data$noise_mag == unqnoise[ii] & data$algorithm%in%unqalg[jj]
      set0 = data$cluster[I0]
      set1 = data$cluster[I1]
      
      # Ai
      tmp = numeric(length(set1))
      for (kk in 1:length(tmp)) {
        cc = cc+1
        ai$algorithm[cc] = unqalg[jj]
        ai$noise_mag[cc] = unqnoise[ii]
        ai$ai[cc] = calcA(set1[kk], set0)
      }
    }
  }
  ai = ai[1:cc,]
  ai$noise_mag = as.numeric(ai$noise_mag)
  
  save("ai", file = sf)
}

I = ai$noise_mag %in% 0 | ai$noise_mag %in% 0.01 | ai$noise_mag %in% 0.02
ggplot(ai[I,], aes(factor(noise_mag), y=ai)) + 
  geom_violin() + geom_jitter(width=.02,alpha=.05) + facet_grid(~algorithm) +
  ylab("Complex reproducibility, A") + xlab("Interactome FDR") + 
  coord_cartesian(ylim=c(0,1)) + theme_bw()
fn = "/Users/Mercy/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2C_v02.pdf"
ggsave(fn,width=10, height=3)
