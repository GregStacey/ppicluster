# Same as full_netw_analysis.R, but analyzes chem-protein and 
# facebook networks (instead of CORUM).
#
######################################################################

source("functions.R")

# read chem-protein and facebook networks
fns = c("../data/ChCh-Miner_durgbank-chem-chem.tsv", "../data/email-Eu-core.txt")
seps = c("\t", " ")
ints = list()
for (ii in 1:length(fns)) {
  ints[[ii]] = as.data.frame(read_delim(fns[ii], delim=seps[ii], quote = ""))
  colnames(ints[[ii]]) = c("protA","protB")
  
  # remove self-interactions
  I = ints[[ii]][,1] == ints[[ii]][,2]
  ints[[ii]] = ints[[ii]][!I,]
}
ints[[2]] = ints[[2]] + 1


# 
noise.range = c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.00)

# load Matlab results
# 1 iterations, MCL, CO, CO+MCL
fn = "../data/clusters_chem_email.txt"
clusters0 = as.data.frame(read_tsv(fn))

# make dummy dataframe to fill
clusters2 = data.frame(network = numeric(10^6),
                       noise_mag = numeric(10^6),
                       algorithm = character(10^6),
                       cluster = character(10^6), stringsAsFactors = F)

# add 1 iteration of PAM + walktrap
cc = 0
networks = c("chem", "email")
n.clusters = c(193, 200)
for (uu in 1:length(fns)) {
  for (ii in 1:length(noise.range)) {
    print(paste("clustering",networks[uu],"at noise=", noise.range[ii]))

    # get shuffled corum
    ints.shuffle = shufflecorum(ints[[uu]], noise.range[ii])
    
    # pam
    print("pam")
    pam.cluster = pamclust(ints.shuffle, n.clusters[uu])
    for (jj in 1:length(pam.cluster)) {
      cc = cc+1
      clusters2$network[cc] = networks[uu]
      clusters2$noise_mag[cc] = noise.range[ii]
      clusters2$algorithm[cc] = "pam"
      clusters2$cluster[cc] = pam.cluster[jj]
    }
    
    # walktrap
    print("walk")
    graph.object = graph_from_edgelist(as.matrix(ints.shuffle), directed = F)
    walk.cluster = walktrap.community(graph.object, steps = 1.5, modularity=T)
    for (jj in 1:length(walk.cluster)) {
      if (length(walk.cluster[[jj]]) < 3) next
      cc = cc+1
      clusters2$network[cc] = networks[uu]
      clusters2$noise_mag[cc] = noise.range[ii]
      clusters2$algorithm[cc] = "walk"
      clusters2$cluster[cc] = paste(walk.cluster[[jj]], collapse=";")
    }
    
    # write in case of crash
    fn = "../data/clusters_facebook_netw_pamwalk.txt"
    write_tsv(clusters2[1:cc,], path=fn)
  }
}
clusters2 = clusters2[1:cc,]

# combine everything
clusters = rbind(clusters0, clusters2)
unqmags = unique(clusters$noise_mag)
unqalgs = unique(clusters$algorithm)
unqnets = unique(clusters$network)


# write in case of crash
fn = "../data/clusters_facebook_netw_pamwalk.txt"
write_tsv(clusters, path=fn)


# calculate Ji1 all clusters (Compare each cluster to its unnoised version)
clusters$Ji1 = numeric(nrow(clusters))
for (ii in 1:length(unqnets)) {
  for (jj in 1:length(unqalgs)) {
    print(paste("algorithm",unqalgs[jj]))
    I0 = clusters$algorithm==unqalgs[jj] & clusters$network==unqnets[ii]
    ref.clusters = clusters$cluster[I0 & clusters$noise_mag==0]
    for (kk in 1:length(unqmags)) {
      print(paste("   noise",unqmags[kk]))
      I = which(I0 & clusters$noise_mag==unqmags[kk])
      these.clusters = clusters$cluster[I]
      for (mm in 1:length(I)) {
        clusters$Ji1[I[mm]] = calcA(these.clusters[mm], ref.clusters)
      }
    }
  }
}



# write
clusters$iter = 1
fn = "../data/clusters_facebook_netw_pamwalk.txt"
write_tsv(clusters, path=fn)


