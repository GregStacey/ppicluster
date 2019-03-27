# Do full network-noise analysis.
# That is:
#   1. Shuffle CORUM x10 iterations
#   2. Cluster iter-shuffled-CORUM with all algorithms
#   3. Compare each cluster to its unnoised version (Ji1)
#   4. Compare each cluster to its iter=1 version (Ji2)
#
# Make this dataframe:
# df = data.frame(cluster = character(10^6),
#                 iter = numeric(10^6),
#                 algorithm = character(10^6),
#                 noise_mag = numeric(10^6),
#                 Ji1 = numeric(10^6),
#                 Ji2 = numeric(10^6), stringsAsFactors = F)
#
# The catch: You don't have three of the algorithms implemented in R!
# Fortunately these were already done in Matlab.
# Therefore 
#   - start from the Matlab results
#   - calculate Ji1 and Ji2
#   - continue 
#
######################################################################

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

# load Matlab results
# 10 iterations, MCL, CO, CO+MCL
fn = "../data/clusters_Ai_vs_fdr.txt"
clusters0 = as.data.frame(read_tsv(fn))

# make dummy dataframe to fill
clusters2 = data.frame(iter = numeric(10^6),
                       noise_mag = numeric(10^6),
                       algorithm = character(10^6),
                       cluster = character(10^6), stringsAsFactors = F)

# add 10 iterations of PAM and hierarchical
iterMax = 10
cc = 0
for (iter in 1:iterMax) {
  for (ii in 1:length(noise.range)) {
    # get shuffled corum
    ints.shuffle = shufflecorum(ints.corum, noise.range[ii])
    
    # hierarchical
    hi.cluster = hiclust(ints.shuffle, 1500)
    for (jj in 1:length(hi.cluster)) {
      cc = cc+1
      clusters2$iter[cc] = iter
      clusters2$noise_mag[cc] = noise.range[ii]
      clusters2$algorithm[cc] = "hierarchical"
      clusters2$cluster[cc] = hi.cluster[jj]
    }
    
    # pam
    pam.cluster = pamclust(ints.shuffle, 1500)
    for (jj in 1:length(hi.cluster)) {
     cc = cc+1
     clusters2$iter[cc] = iter
     clusters2$noise_mag[cc] = noise.range[ii]
     clusters2$algorithm[cc] = "pam"
     clusters2$cluster[cc] = hi.cluster[jj]
    }
  }
}
clusters2 = clusters2[1:cc,]

# combine everything
clusters = rbind(clusters0, clusters2)
unqiters = unique(clusters$iter)
unqmags = unique(clusters$noise_mag)
unqalgs = unique(clusters$algorithm)


# calculate Ji1 all clusters (Compare each cluster to its unnoised version)
clusters$Ji1 = numeric(nrow(clusters))
for (ii in 1:length(unqiters)) {
  print(paste("Ji1: iter",ii))
  for (jj in 1:length(unqalgs)) {
    print(paste("      algorithm",unqalgs[jj]))
    I0 = clusters$iter==unqiters[ii] & clusters$algorithm==unqalgs[jj]
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
fn = "../data/clusters_full_netw.txt"
write_tsv(clusters, path=fn)

# calculate Ji2 all clusters (Compare each cluster to its iter=1 version)
clusters$Ji2 = numeric(nrow(clusters))
for (kk in 1:length(unqmags)) {
  print(paste("Ji2: noise",unqmags[kk]))
  for (jj in 1:length(unqalgs)) {
    print(paste("      algorithm",unqalgs[jj]))
    I0 = clusters$algorithm==unqalgs[jj] & clusters$noise_mag==unqmags[kk]
    ref.clusters = clusters$cluster[I0 & clusters$iter==1]
    for (ii in 1:length(unqiters)) {
      print(paste("        iter",unqmags[kk]))
      I = which(I0 & clusters$iter==unqiters[ii])
      these.clusters = clusters$cluster[I]
      for (mm in 1:length(I)) {
        clusters$Ji2[I[mm]] = calcA(these.clusters[mm], ref.clusters)
      }
    }
  }
}
# write
fn = "../data/clusters_full_netw.txt"
write_tsv(clusters, path=fn)


