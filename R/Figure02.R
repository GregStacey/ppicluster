
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
for (uu in c(8)) {
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
  
  # plot
  ggplot() + geom_scatterpie(aes(x=x, y=y, r=sqrt(size)), data=df, cols=names(df)[1:3]) +
    #scale_fill_manual(values=c("#4dac26","#d01c8b","#dddddd")) + # (green, red, grey) = (grey, black, white)
    scale_fill_manual(values=c("grey", "white","black")) + # (green, red, grey) = (grey, black, white)
    coord_fixed() + theme_void() + theme(legend.position="none")
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_", uu, ".pdf", sep="")
  ggsave(sf, width=10, height=5)
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_", uu, ".png", sep="")
  ggsave(sf, width=10, height=5, dpi=1000)
  
  # zoom in
  ggplot() + geom_scatterpie(aes(x=x, y=y, r=sqrt(size)), data=df, cols=names(df)[1:3]) + 
    #scale_fill_manual(values=c("#4dac26","#d01c8b","#dddddd")) +
    scale_fill_manual(values=c("grey", "white","black")) + # (green, red, grey) = (grey, black, white)
    coord_fixed() + theme_void() + theme(legend.position="none") + coord_cartesian(xlim = c(80, 130))
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_zoom_", uu, ".pdf", sep="")
  ggsave(sf, width=2, height = 5)
  sf = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/figure02/pies_zoom_", uu, ".png", sep="")
  ggsave(sf, width=2, height = 5, dpi=1000)
}



sf = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/data/fig2B_v06.Rda"
if (T){
  load(sf)
} else {

  # corum
  fn = "../data/clusters_full_netw_walktrap2.txt"
  Ji2 = as.data.frame(read_tsv(fn))
  Ji2 = Ji2[Ji2$iter==1,]
  Ji2$network = "corum"
  Ji2 = Ji2[,c("iter", "network", "noise_mag", "algorithm", "cluster", "Ji1")]
  
  # email+chem
  fn = "../data/clusters_facebook_netw_pamwalk.txt"
  Ji = as.data.frame(read_tsv(fn))
  Ji = Ji[,c("iter", "network", "noise_mag", "algorithm", "cluster", "Ji1")]
  
  # merge
  Ji = rbind(Ji, Ji2)
  
  Ji$network[Ji$network=="chem"] = "DrugBank"
  Ji$network[Ji$network=="corum"] = "CORUM"
  Ji$network[Ji$network=="email"] = "email-Eu"
  
  # make cluster.size, remove all clusters with size<3
  Ji$cluster.size = sapply((sapply(Ji$cluster, strsplit, ";")), length)
  Ji = Ji[Ji$cluster.size>2,]
  
  I = !Ji$algorithm == "hierarchical" & Ji$iter==1
  Ji = Ji[I,]
  unqnets = unique(Ji$network)
  unqalg = unique(Ji$algorithm)
  unqnoise = unique(Ji$noise_mag)
  sim = data.frame(network = character(10^4),
                   measure = character(10^4),
                   algorithm = character(10^4),
                   noise_mag = character(10^4),
                   Ji1 = numeric(10^4), stringsAsFactors = F)
  cc = 0
  cc2 = 0
  for (uu in 1:length(unqnets)) {
    for (ii in 1:length(unqnoise)) {
      for (jj in 1:length(unqalg)) {
      
        I1 = Ji$noise_mag == unqnoise[ii] & Ji$algorithm%in%unqalg[jj] & Ji$network==unqnets[uu]
        if (sum(I1)==0) next
        
        # Ji1
        cc = cc+1
        sim$measure[cc] = "ai"
        sim$network[cc] = unqnets[uu]
        sim$algorithm[cc] = unqalg[jj]
        sim$noise_mag[cc] = unqnoise[ii]
        sim$Ji1[cc] = mean(Ji$Ji1[I1])
      }
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
Ji$algorithm[Ji$algorithm=="walk"] = "walktrap"
sim$algorithm[sim$algorithm=="co_mcl"] = "CO+MCL"
sim$algorithm[sim$algorithm=="co"] = "CO"
sim$algorithm[sim$algorithm=="mcl"] = "MCL"
sim$algorithm[sim$algorithm=="pam"] = "k-Med"
sim$algorithm[sim$algorithm=="walk"] = "walktrap"
sim$Ji1 = sim$Ji




# 2B ##### ------------------------------------------- #####
# line plots of J_i vs noise

I = sim$measure=="ai" & sim$measure=="ai" & sim$network=="CORUM"
ggplot(sim[I,], aes(x=noise_mag, y=Ji1)) + geom_line() + 
  facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
  geom_point(data = Ji[Ji$network=="CORUM",], alpha=0.025, color="black") +
  ylab("Similarity to un-noised clusters (Ji)") + xlab("Interactome FPR")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2B_v08.png"
ggsave(fn,width=10, height=2.6)



# 2C ##### ------------------------------------------- #####
# violin plots of Ji1 vs noise

I = (Ji$noise_mag %in% 0 | Ji$noise_mag %in% 0.01 | Ji$noise_mag %in% 0.02) & Ji$network=="CORUM"
ggplot(Ji[I,], aes(factor(noise_mag), y=Ji1)) + 
  geom_violin() + geom_jitter(width=.02,alpha=.05) + facet_grid(~algorithm) +
  ylab("Similarity to un-noised clusters (Ji)") + xlab("Interactome FPR") + 
  coord_cartesian(ylim=c(0,1)) + theme_bw()
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2C_v06.pdf"
ggsave(fn,width=10, height=2.6)
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2C_v06.png"
ggsave(fn,width=10, height=2.6)


# 2D ##### ------------------------------------------- #####
# noise-amplification

I = sim$measure=="ai"& sim$network=="CORUM"
sim$noise.amplification = (1-sim$Ji1) / sim$noise_mag
ggplot(sim[I,], aes(x=noise_mag, y=(noise.amplification))) + geom_line() + 
  facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
  geom_hline(yintercept = 1, linetype="dashed") +
  ylab("Error amplification by clustering") + xlab("Interactome FPR") +
  scale_y_continuous(breaks = c(1,10,20),
                     labels=c("1x","10x", "20x"))
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2D_v03.pdf"
ggsave(fn,width=10, height=2.6)


# 2E,F ##### ------------------------------------------- #####

I = sim$measure=="ai" & sim$measure=="ai" & sim$network %in% c("DrugBank","email-Eu") & sim$algorithm=="MCL"
ggplot(sim[I,], aes(x=noise_mag, y=Ji1)) + geom_line() + 
  theme_bw() + theme(legend.position="none") + facet_grid(~network) +
  geom_point(data = Ji[Ji$network%in% c("DrugBank","email-Eu") &Ji$algorithm=="MCL",], alpha=0.1, color="black") +
  ylab("Similarity to un-noised clusters (Ji)") + xlab("Interactome FPR")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2E_v01.png"
ggsave(fn,width=4.2, height=2.6)

I = (Ji$noise_mag %in% 0 | Ji$noise_mag %in% 0.01 | Ji$noise_mag %in% 0.02) & Ji$network%in% c("DrugBank","email-Eu") & 
  Ji$algorithm=="MCL"
ggplot(Ji[I,], aes(factor(noise_mag), y=Ji1)) + 
  geom_violin() + geom_jitter(width=.02,alpha=.1) + facet_grid(~network) +
  ylab("Similarity to un-noised (Ji)") + xlab("Interactome FPR") + 
  coord_cartesian(ylim=c(0,1)) + theme_bw()
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2F_v01.pdf"
ggsave(fn,width=4.2, height=2.6)
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_2F_v01.png"
ggsave(fn,width=4.2, height=2.6)




# supp of drugbank and email sets ##### ------------------------------------------- #####
for (ss in 1:2) {
  sets = c("DrugBank","email-Eu")
  # line plots of J_i vs noise
  I = sim$measure=="ai" & sim$measure=="ai" & sim$network==sets[ss]
  ggplot(sim[I,], aes(x=noise_mag, y=Ji1)) + geom_line() + 
    facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
    geom_point(data = Ji[Ji$network==sets[ss],], alpha=0.025, color="black") +
    ylab("Similarity to un-noised clusters (Ji)") + xlab("Interactome FPR")
  fn = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1A",ss,".png", sep="")
  ggsave(fn,width=10, height=2.6)
  
  # violin plots of Ji1 vs noise
  I = (Ji$noise_mag %in% 0 | Ji$noise_mag %in% 0.01 | Ji$noise_mag %in% 0.02) & Ji$network==sets[ss]
  ggplot(Ji[I,], aes(factor(noise_mag), y=Ji1)) + 
    geom_violin() + geom_jitter(width=.02,alpha=.05) + facet_grid(~algorithm) +
    ylab("Similarity to un-noised clusters (Ji)") + xlab("Interactome FPR") + 
    coord_cartesian(ylim=c(0,1)) + theme_bw()
  fn = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1B",ss,".png", sep="")
  ggsave(fn,width=10, height=2.6)
  fn = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1B",ss,".pdf", sep="")
  ggsave(fn,width=10, height=2.6)
  
  # noise-amplification
  I = sim$measure=="ai"& sim$network==sets[ss]
  sim$noise.amplification = (1-sim$Ji1) / sim$noise_mag
  ggplot(sim[I,], aes(x=noise_mag, y=(noise.amplification))) + geom_line() + 
    facet_grid(~algorithm) + theme_bw() + theme(legend.position="none") + 
    geom_hline(yintercept = 1, linetype="dashed") +
    ylab("Error amplification by clustering") + xlab("Interactome FPR") +
    scale_y_continuous(breaks = c(1,10,20),
                       labels=c("1x","10x", "20x"))
  fn = paste("/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup1C",ss,".pdf", sep="")
  ggsave(fn,width=10, height=2.6)
}

