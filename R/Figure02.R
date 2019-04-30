
require(scatterpie)
source("functions.R")


fnsave = "../data/clusters_wshuffle_coexp.Rda"
load(fnsave)
clusters$noise_mag = as.numeric(clusters$noise_mag)

fnsave = "../data/cluster3.txt"
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
  pci = matrix(nrow=length(set1), ncol=3)#length(set0))
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
    nn.notintersect.set1 = length(cluster1) - nn.intersect
    nn.notintersect.set0 = length(cluster0) - nn.intersect
    
    pci[ii,] = c(nn.intersect, nn.notintersect.set1, nn.notintersect.set0)
  }
  pci = as.data.frame(pci)
  pci$cluster = 1:nrow(pci)
  pci$size = unlist(lapply(lapply(lapply(set1, strsplit, ";"), "[[", 1), length))

  # tile x and y
  # assume that cluster with sqrt(size) needs a box of width=m*sqrt(size)
  gridx = 250
  df = pci[,]
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
  df = df[,c(1:3,ncol(df) + seq(from=-3, to=0, by=1))]
  df[df==0] = 10^-4
  
  plot
  ggplot() + geom_scatterpie(aes(x=x, y=y, r=sqrt(size)), data=df, cols=names(df)[1:3]) +
   scale_fill_manual(values=c("#4dac26","#d01c8b","#dddddd")) +
   coord_fixed() + theme_void() + theme(legend.position="none")
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_", uu, ".pdf", sep="")
  ggsave(sf, width=10, height=5)
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_", uu, ".png", sep="")
  ggsave(sf, width=10, height=5, dpi=1000)
  
  # zoom in
  ggplot() + geom_scatterpie(aes(x=x, y=y, r=sqrt(size)), data=df, cols=names(df)[1:3]) + 
    scale_fill_manual(values=c("#4dac26","#d01c8b","#dddddd")) +
    coord_fixed() + theme_void() + theme(legend.position="none") + coord_cartesian(xlim = c(80, 130))
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_zoom_", uu, ".pdf", sep="")
  ggsave(sf, width=2, height = 5)
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_zoom_", uu, ".png", sep="")
  ggsave(sf, width=2, height = 5, dpi=1000)
}



sf = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/data/fig2B_v04.Rda"
if (F){
  load(sf)
} else {
  fn = "../data/clusters_full_netw.txt"
  Ji = as.data.frame(read_tsv(fn))
  
  # make cluster.size, remove all clusters with size<3
  Ji$cluster.size = sapply((sapply(Ji$cluster, strsplit, ";")), length)
  Ji = Ji[Ji$cluster.size>2,]
  
  I = !Ji$algorithm == "hierarchical" & Ji$iter==1
  Ji = Ji[I,]
  unqalg = unique(Ji$algorithm)
  unqnoise = unique(Ji$noise_mag)
  sim = data.frame(measure = character(10^4),
                   algorithm = character(10^4),
                   noise_mag = character(10^4),
                   Ji1 = numeric(10^4), stringsAsFactors = F)
  cc = 0
  cc2 = 0
  for (ii in 1:length(unqnoise)) {
    print(unqnoise[ii])
    for (jj in 1:length(unqalg)) {
      print(paste("  ",unqalg[jj]))
      I0 = Ji$noise_mag == 0 & Ji$algorithm%in%unqalg[jj]
      I1 = Ji$noise_mag == unqnoise[ii] & Ji$algorithm%in%unqalg[jj]
      set0 = Ji$cluster[I0]
      set1 = Ji$cluster[I1]
      if (sum(I1)==0 | sum(I0)==0) next
      
      # geomacc
      cc = cc+1
      sim$measure[cc] = "ga"
      sim$algorithm[cc] = unqalg[jj]
      sim$noise_mag[cc] = unqnoise[ii]
      sim$Ji1[cc] = geomacc(set0, set1)
      
      # matching ratio
      cc = cc+1
      sim$measure[cc] = "mmr"
      sim$algorithm[cc] = unqalg[jj]
      sim$noise_mag[cc] = unqnoise[ii]
      sim$Ji1[cc] = matchingratio(set0, set1)
      
      # nmi
      cc = cc+1
      unqprots = c(unlist(lapply(set0,strsplit,";")),unlist(lapply(set1,strsplit,";")))
      # prepare set0
      df0 = data.frame(X=1:length(unqprots),Y=numeric(length(unqprots)))
      for (kk in 1:length(unqprots)) {
        ia = grep(unqprots[kk], set0)
        if (length(ia)==0) next
        df0$Y[kk] = ia[1]
      }
      df0 = df0[!df0$Y==0,]
      # prepare set1
      df1 = data.frame(X=1:length(unqprots),Y=numeric(length(unqprots)))
      for (kk in 1:length(unqprots)) {
        ia = grep(unqprots[kk], set1)
        if (length(ia)==0) next
        df1$Y[kk] = ia[1]
      }
      df1 = df1[!df1$Y==0,]
      sim$measure[cc] = "NMI"
      sim$algorithm[cc] = unqalg[jj]
      sim$noise_mag[cc] = unqnoise[ii]
      sim$Ji1[cc] = NMI(df0, df1)[[1]]
      
      
      # Ji1
      cc = cc+1
      sim$measure[cc] = "ai"
      sim$algorithm[cc] = unqalg[jj]
      sim$noise_mag[cc] = unqnoise[ii]
      # tmp = numeric(length(set1))
      # for (kk in 1:length(tmp)) {
      #   cc2 = cc2+1
      #   tmp[kk] = calcA(set1[kk], set0)
      #   Ji$ji[cc2] = tmp[kk]
      #   Ji$algorithm[cc2] = unqalg[jj]
      #   Ji$noise_mag[cc2] = unqnoise[ii]
      #   Ji$cluster.size[cc2] = length(unlist(strsplit(set1[kk], ";")))
      # }
      sim$Ji1[cc] = mean(Ji$Ji1[I1])
    }
  }
  sim = sim[1:cc,]
  sim$noise_mag = as.numeric(sim$noise_mag)
  sim$Ji[sim$Ji==0] = NA
  Ji$noise_mag = as.numeric(Ji$noise_mag)

  save(list=c("sim","Ji"), file = sf)
}

Ji$measure = Ji$algorithm
Ji$algorithm[Ji$algorithm=="co_mcl"] = "CO+MCL"
Ji$algorithm[Ji$algorithm=="co"] = "CO"
Ji$algorithm[Ji$algorithm=="mcl"] = "MCL"
Ji$algorithm[Ji$algorithm=="pam"] = "k-Med"
sim$algorithm[sim$algorithm=="co_mcl"] = "CO+MCL"
sim$algorithm[sim$algorithm=="co"] = "CO"
sim$algorithm[sim$algorithm=="mcl"] = "MCL"
sim$algorithm[sim$algorithm=="pam"] = "k-Med"
sim$Ji1 = sim$Ji





# 2B ##### ------------------------------------------- #####
# line plots of J_i vs noise

I = sim$measure=="ai"
ggplot(sim[I,], aes(x=noise_mag, y=Ji1)) + geom_line() + 
  facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
  geom_point(data = Ji, alpha=0.025, color="black") +
  ylab("Similarity to un-noised clusters (Ji)") + xlab("Interactome FPR")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2B_v05.pdf"
ggsave(fn,width=10, height=3)
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2B_v05.png"
ggsave(fn,width=10, height=3)



# 2C ##### ------------------------------------------- #####
# violin plots of Ji1 vs noise

I = Ji$noise_mag %in% 0 | Ji$noise_mag %in% 0.01 | Ji$noise_mag %in% 0.02
ggplot(Ji[I,], aes(factor(noise_mag), y=Ji1)) + 
  geom_violin() + geom_jitter(width=.02,alpha=.05) + facet_grid(~algorithm) +
  ylab("Similarity to un-noised (Ji)") + xlab("Interactome FPR") + 
  coord_cartesian(ylim=c(0,1)) + theme_bw()
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2C_v03.pdf"
ggsave(fn,width=10, height=3)
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2C_v03.png"
ggsave(fn,width=10, height=3)


# 2D ##### ------------------------------------------- #####
# noise-amplification

I = sim$measure=="ai"
sim$noise.amplification = (1-sim$Ji1) / sim$noise_mag
ggplot(sim[I,], aes(x=noise_mag, y=(noise.amplification))) + geom_line() + 
  facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
  geom_hline(yintercept = 1, linetype="dashed") +
  ylab("Error amplification by clustering") + xlab("Interactome FPR") +
  scale_y_continuous(breaks = c(1,10,20),
                     labels=c("1x","10x", "20x"))
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2D_v01.pdf"
ggsave(fn,width=10, height=3)




cn = data.frame(year = 2008:2013,
                US = c(799, 724,677,593,544,541),
                Canada = c(157, 164, 199, 206, 197, 197),
                International = c(186, 184, 242, 206, 146, 145))
df = melt(cn, id.vars = "year")
ggplot(df, aes(x=year, y=value, color=variable)) + geom_line() +
  xlab("Year") + ylab("Number of\nCoffee News franchises") + ggtitle("Why did I make this")
ggsave("/Users/gregstacey/Desktop//coffee_news_trends.png",
       width=5, height=3)


cn.pc = data.frame(year = 2008:2013,
                US = c(799, 724,677,593,544,541) / 346,
                Canada = c(157, 164, 199, 206, 197, 197) / 37)
df = melt(cn.pc, id.vars = "year")
ggplot(df, aes(x=year, y=value, color=variable)) + geom_line() +
  xlab("Year") + ylab("Number of\nCoffee News franchises,\nper million pop.") + 
  ggtitle("Okay this is cool") + ylim(c(0,6))
ggsave("/Users/gregstacey/Desktop//coffee_news_trends.png",
       width=5, height=3)


