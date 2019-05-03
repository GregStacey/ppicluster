
source("unctions.R")

fn = "../data/clusters_full_netw.txt"
Ji = as.data.frame(read_tsv(fn))

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
                  sum(Ji$Ji2[I]>0.65, na.rm=T) / sum(!is.na(Ji$Ji2[I])),
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
  xlab("Interactome FPR") + ylab("Fraction of clusters\nreproducible (Ai>0.65)") + 
  coord_cartesian(ylim=c(0,1)) + theme_bw() 
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_5A_v01.png"
ggsave(fn,width=10, height=3)



I = Ji$iter>1 & Ji$algorithm=="CO"
a0 = summary(glm(formula = Ji2 ~ noise_mag  + cluster.size, data = Ji[I,]))
a = summary(glm(formula = Ji2 ~ noise_mag, data = Ji[I,]))
c("CO",1 - a$deviance/a$null.deviance,1 - a0$deviance/a0$null.deviance)

I = Ji$iter>1 & Ji$algorithm=="CO+MCL" 
a0 = summary(glm(formula = Ji2 ~ noise_mag  + cluster.size, data = Ji[I,]))
a = summary(glm(formula = Ji2 ~ noise_mag, data = Ji[I,]))
c("CO+MCL",1 - a$deviance/a$null.deviance,1 - a0$deviance/a0$null.deviance)

I = Ji$iter>1 & Ji$algorithm=="k-Med" 
a0 = summary(glm(formula = Ji2 ~ noise_mag  + cluster.size, data = Ji[I,]))
a = summary(glm(formula = Ji2 ~ noise_mag, data = Ji[I,]))
c("k-Med",1 - a$deviance/a$null.deviance,1 - a0$deviance/a0$null.deviance)

I = Ji$iter>1 & Ji$algorithm=="MCL" 
a0 = summary(glm(formula = Ji2 ~ noise_mag  + cluster.size, data = Ji[I,]))
a = summary(glm(formula = Ji2 ~ noise_mag, data = Ji[I,]))
c("MCL",1 - a$deviance/a$null.deviance,1 - a0$deviance/a0$null.deviance)


