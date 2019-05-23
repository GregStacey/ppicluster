
source("functions.R")

fn = "../data/clusters_full_netw.txt"
Ji = as.data.frame(read_tsv(fn))
Ji = Ji[!Ji$algorithm=="hierarchical",]
Ji$measure = Ji$algorithm
Ji$algorithm[Ji$algorithm=="co_mcl"] = "CO+MCL"
Ji$algorithm[Ji$algorithm=="co"] = "CO"
Ji$algorithm[Ji$algorithm=="mcl"] = "MCL"
Ji$algorithm[Ji$algorithm=="pam"] = "k-Med"

# make cluster.size, remove all clusters with size<3
Ji$cluster.size = sapply((sapply(Ji$cluster, strsplit, ";")), length)
Ji = Ji[Ji$cluster.size>2,]

# make size.factor
Ji$noise_mag = as.numeric(Ji$noise_mag)
Ji$size.factor = character(nrow(Ji))
Ji$size.factor[Ji$cluster.size<=3] = "N<=3"
Ji$size.factor[Ji$cluster.size<=6 & Ji$cluster.size>3] = "3<N<=6"
Ji$size.factor[Ji$cluster.size<=12 & Ji$cluster.size>6] = "6<N<=12"
Ji$size.factor[Ji$cluster.size>12 & Ji$cluster.size<50] = "N>12"
Ji$size.factor[Ji$cluster.size>=50] = "Big"
Ji$size.factor = factor(Ji$size.factor, levels= c("N<=3", "3<N<=6", "6<N<=12", "N>12", "Big"))

# 5A. % reproducible vs FPR
nn = 10^3
dm = data.frame(x = numeric(nn), 
                y = numeric(nn),
                y2 = numeric(nn),
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
      cc = cc+1
      I = Ji$noise_mag==this.noise & 
        Ji$algorithm%in%this.algorithm & 
        Ji$size.factor%in%this.size & Ji$iter>1
      dm[cc,] = c(this.noise, mean(Ji$Ji2[I], na.rm=T), 
                  sum(Ji$Ji2[I]>0.5, na.rm=T) / sum(!is.na(Ji$Ji2[I])),
                  this.algorithm, this.size)
    }
  }
}
dm = dm[1:cc,]
dm$x = as.numeric(dm$x)
dm$y = as.numeric(dm$y)
dm$y2 = as.numeric(dm$y2)
tmp = sort(unique(Ji$size.factor))
dm$size.factor = tmp[as.numeric(dm$size.factor)]
dm = dm[!is.na(dm$y2),]


ggplot(dm, aes(x=x, y=y2, color=size.factor)) + 
  geom_line(alpha=0.6, size=2) +  facet_grid(~algorithm) +
  xlab("Interactome FPR") + ylab("Fraction of clusters\nreproducible (Ji>0.5)") + 
  coord_cartesian(ylim=c(0,1)) + theme_bw() + scale_color_grey() + 
  theme(legend.position = "none")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_5A_v01.png"
ggsave(fn,width=10, height=3)


