# add supplementary results - add+remove noise is also amplified by clustering


source("functions.R")

# binarize corum
fn = "../data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))
corum = corum[corum$Organism=="Human",]

# turn corum into a list of protein pairs
ints.corum = data.frame(protA = character(10^6),
                        protB = character(10^6), stringsAsFactors = F)
cc = 0
for (ii in 1:nrow(corum)) {
  print(ii)
  prots = sort(unlist(strsplit(corum$`subunits(UniProt IDs)`[ii], ";")))
  if (length(prots)<2) next
  pairs.prots = t(combn(prots, 2))
  
  I = (cc+1) : (cc+nrow(pairs.prots))
  ints.corum$protA[I] = pairs.prots[,1]
  ints.corum$protB[I] = pairs.prots[,2]
  cc = cc+length(I)
}
ints.corum = ints.corum[1:cc,]
ints.corum = distinct(ints.corum)

# 
noise.range = c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.00)



# load (shuffle+add+remove) * (three algorithms)
fn = "../data/clusters_coexp.Rda" # clusters
load(fn)
clusters = clusters[clusters$data_type=="corum", 
                    c("data_type", "noise_type", "noise_mag", "algorithm", "cluster")]
names(clusters) = c("network", "noise_type", "noise_mag", "algorithm", "cluster")


# make dummy dataframe to fill
clusters2 = data.frame(network = rep("corum", 10^6),
                       noise_mag = numeric(10^6),
                       noise_type = character(10^6),
                       algorithm = character(10^6),
                       cluster = character(10^6), stringsAsFactors = F)

# add 1 iterations of PAM and walktrap
cc = 0
noise.types = c("network_remove","network_add")
for (hh in 1:length(noise.types)) {
  for (ii in 1:length(noise.range)) {
    print(paste("clustering",noise.types[hh],"at noise =", noise.range[ii]))
    
    # get shuffled corum
    if (noise.types[hh] == "network_add") ints.noised = addcorum(ints.corum, noise.range[ii])
    if (noise.types[hh] == "network_remove") ints.noised = removecorum(ints.corum, noise.range[ii])
    
    # if you removed 100%... change it to 95%
    if (nrow(ints.noised)==0 & noise.range[ii] & noise.types[hh]=="network_remove"){
      ints.noised = removecorum(ints.corum, .95)
    }
    
    # walktrap
    graph.object = graph_from_edgelist(as.matrix(ints.noised), directed = F)
    walk.cluster = walktrap.community(graph.object)
    for (jj in 1:length(walk.cluster)) {
      if (length(walk.cluster[[jj]]) < 3) next
      cc = cc+1
      clusters2$noise_type[cc] = noise.types[hh]
      clusters2$noise_mag[cc] = noise.range[ii]
      clusters2$algorithm[cc] = "walk"
      clusters2$cluster[cc] = paste(walk.cluster[[jj]], collapse=";")
    }
    
    # pam
    pam.cluster = pamclust(ints.noised, 1500)
    for (jj in 1:length(pam.cluster)) {
      cc = cc+1
      clusters2$noise_type[cc] = noise.types[hh]
      clusters2$noise_mag[cc] = noise.range[ii]
      clusters2$algorithm[cc] = "pam"
      clusters2$cluster[cc] = pam.cluster[jj]
    }
  }
}
clusters2 = clusters2[1:cc,]

# combine everything
clusters = rbind(clusters, clusters2)
unqnoise = unique(clusters$noise_type)
unqmags = unique(clusters$noise_mag)
unqalgs = unique(clusters$algorithm)


# calculate Ji1 all clusters (Compare each cluster to its unnoised version)
clusters$Ji1 = rep(NA, nrow(clusters))
for (ii in 1:length(unqnoise)) {
  print(paste("Ji1: noise",unqnoise[ii]))
  for (jj in 1:length(unqalgs)) {
    print(paste("      algorithm",unqalgs[jj]))
    I0 = clusters$noise_type==unqnoise[ii] & clusters$algorithm==unqalgs[jj]
    ref.clusters = clusters$cluster[I0 & clusters$noise_mag==0]
    for (kk in 1:length(unqmags)) {
      print(paste("        noise",unqmags[kk]))
      I = which(I0 & clusters$noise_mag==unqmags[kk])
      these.clusters = clusters$cluster[I]
      for (mm in 1:length(I)) {
        clusters$Ji1[I[mm]] = calcA(these.clusters[mm], ref.clusters)
      }
    }
  }
}
# write
fn = "../data/clusters_full_netw_walktrap_addremove.txt"
write_tsv(clusters, path=fn)

