
source("functions.R")


# matlab results (chromatograms + network, three algorithms)
fnsave = "../data/clusters_wshuffle_coexp.Rda"
load(fnsave)
clusters$noise_mag = as.numeric(clusters$noise_mag)

# pam
fnsave = "../data/cluster3.txt"
clusters2 = as.data.frame(read_tsv(fnsave))
clusters2$noise_mag = as.numeric(clusters2$noise_mag)
clusters2 = clusters2[clusters2$data_type=="corum",]

# walktrap
fn = "../data/clusters_full_netw_walktrap2.txt"
clusters3 = as.data.frame(read_tsv(fn))
clusters3 = clusters3[clusters3$iter==1,]
clusters3$data_type = "corum"
clusters3$noise_type = "network_shuffle"
clusters3 = clusters3[,c("data_type", "noise_type", "noise_mag", "algorithm", "cluster")]

clusters = rbind(clusters[,1:5], clusters2, clusters3)


# Count the rearrangements. "Simple" analysis
sf = "../data/dfchange.Rda"
if (T) {
  load(sf)
} else {
  df.change = data.frame(noise = numeric(10^4),
                         algorithm = character(10^4),
                         n.new = numeric(10^4), 
                         n.lost = numeric(10^4),
                         n.clusters = numeric(10^4),
                         avg.gained = numeric(10^4),
                         avg.lost = numeric(10^4), stringsAsFactors = F)
  
  unqnoise = unique(data$noise_mag)
  unqalgs = unique(data$algorithm)
  unqalgs = unqals[!unqalgs == "hierarchical"]
  cc = 0
  for (mm in 1:length(unqalgs)) {
    I = clusters$data_type%in%"corum" & clusters$noise_type%in%"network_shuffle"
    data = clusters[I,]
    set0 = data$cluster[data$noise_mag==0 & data$algorithm%in%unqalgs[mm]]
    print(unqalgs[mm])
    for (uu in 1:length(unqnoise)) {
      print(unqnoise[uu])
      set1 = data$cluster[data$noise_mag==unqnoise[uu] & data$algorithm%in%unqalgs[mm]]
      
      # scan set1 for n.gained
      n.gained = numeric(length(set1))
      n.set1 = numeric(length(set1))
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
        n.gained[ii] = sum(!cluster1 %in% cluster0)
        n.set1[ii] = length(cluster1)
      }
      
      # scan set0 for n.lost
      n.lost = numeric(length(set0))
      n.set0 = numeric(length(set0))
      for (ii in 1:length(set0)) {
        # find best-matching cluster
        cluster0 = unlist(strsplit(set0[ii], ";"))
        JJ = numeric(length(set1))
        for (jj in 1:length(set1)) {
          cluster1 = unlist(strsplit(set1[jj], ";"))
          nn.jacc = length(unique(c(cluster1, cluster0)))
          JJ[jj] = length(intersect(cluster1, cluster0)) / nn.jacc 
        }
        I.max = which.max(JJ)[1]
        cluster1 = unlist(strsplit(set1[I.max], ";"))
        n.lost[ii] = sum(!cluster0 %in% cluster1)
        n.set0[ii] = length(cluster0)
      }
      
      cc = cc+1
      df.change$noise[cc] = unqnoise[uu]
      df.change$algorithm[cc] = unqalgs[mm]
      df.change$n.lost[cc] = sum((n.lost-n.set0)==0)
      df.change$n.new[cc] = sum((n.gained-n.set1)==0)
      df.change$n.clusters[cc] = length(set1)
      df.change$avg.gained[cc] = mean(n.gained)
      df.change$avg.lost[cc] = mean(n.lost)
    }
  }
  df.change = df.change[1:cc,]
}
noise.range = c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.00)
df.change = df.change[df.change$noise %in% noise.range, ]



# number of clusters
ggplot(df.change, aes(x=noise, y=n.clusters, color=algorithm)) + geom_line() +
  scale_colour_grey() + xlab("Interactome FPR") + ylab("Number of clusters") +  
  theme_bw()
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup3E1_v02.png"
ggsave(fn,width=4.5, height=3)

# avg gained and lost proteins
df = melt(df.change, id.vars=c("noise", "algorithm"), measure.vars = c("avg.gained", "avg.lost"))
df$variable = as.character(df$variable)
df[df=="avg.gained"] = "Gained"
df[df=="avg.lost"] = "Lost"
ggplot(df, aes(x=noise, y=value, color=algorithm)) + geom_line() + 
  facet_wrap(~variable, scales="free_y") + theme_bw() +
  scale_colour_grey() + xlab("Interactome FPR") + 
  ylab("Number of gained or lost proteins\nper cluster (avg.)") +
  theme(legend.position = "none")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup3E2_v02.png"
ggsave(fn,width=6.5, height=3)

# number of gained and lost clusters
df = melt(df.change, id.vars=c("noise", "algorithm"), measure.vars = c("n.new", "n.lost"))
df$variable = as.character(df$variable)
df[df=="n.new"] = "Gained"
df[df=="n.lost"] = "Lost"
ggplot(df, aes(x=noise, y=value, color=algorithm)) + geom_line() + 
  facet_wrap(~variable, scales="free_y") +
  scale_colour_grey() + xlab("Interactome FPR") + ylab("Number of gained or lost clusters")

