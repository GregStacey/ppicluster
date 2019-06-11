
# B - number and size of clusters

source("functions.R")

fn = "../data/clusters_full_netw_walktrap.txt"
Ji = as.data.frame(read_tsv(fn))
Ji = Ji[!Ji$algorithm=="hierarchical",]
Ji = Ji[Ji$noise_mag %in% c(0,0.01,0.02,0.05,0.1,0.15,0.25,0.5,1),]

Ji$algorithm[Ji$algorithm=="co_mcl"] = "CO+MCL"
Ji$algorithm[Ji$algorithm=="co"] = "CO"
Ji$algorithm[Ji$algorithm=="mcl"] = "MCL"
Ji$algorithm[Ji$algorithm=="pam"] = "k-Med"
Ji$algorithm[Ji$algorithm=="walk"] = "walktrap"

# make cluster.size, remove all clusters with size<3
Ji$cluster.size = sapply((sapply(Ji$cluster, strsplit, ";")), length)
Ji = Ji[Ji$cluster.size>2,]

unqiter = unique(Ji$iter)
unqalgs = c("CO+MCL", "CO", "MCL", "k-Med", "walktrap")
unqmags = c(0,0.01,0.02,0.05,0.1,0.15,0.25,0.5,1)
df = data.frame(noise_mag = numeric(10^3),
                iter = numeric(10^3),
                nsize = numeric(10^3),
                nn = numeric(10^3),
                algorithm = character(10^3), stringsAsFactors = F)
cc = 0
for (kk in 1:length(unqmags)) {
  for (ii in 1:length(unqalgs)) {
    for (jj in 1:length(unqiter)) {
      I = Ji$algorithm==unqalgs[ii] & Ji$iter==unqiter[jj] & Ji$noise_mag==unqmags[kk]

      cc = cc+1
      df$iter[cc] = unqiter[jj]
      df$noise_mag[cc] = unqmags[kk]
      df$nsize[cc] = mean(Ji$cluster.size[I], na.rm=T) # size
      df$nn[cc] = sum(I) # number
      df$algorithm[cc] = unqalgs[ii]
    }
  }
}
df = df[1:cc,]


# size
ggplot(df, aes(x=noise_mag, y=nn)) + geom_point(alpha=.5) + facet_grid(~algorithm)


# number

