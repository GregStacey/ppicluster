# add supplementary results - add AND remove noise is also amplified by clustering
# differerent from full_netw_analysis_addremove.R, because this does add AND remove,
# not add OR remove


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
add.range = c(0, 0.01, 0.05, 0.1, 0.25)
remove.range = c(0, 0.01, 0.05, 0.1, 0.25)


# make dummy dataframe to fill
clusters = data.frame(network = rep("corum", 10^6),
                      add_mag = numeric(10^6),
                      remove_mag = numeric(10^6),
                      algorithm = character(10^6),
                      cluster = character(10^6), stringsAsFactors = F)

# add 1 iterations of PAM and walktrap
cc = 0
for (hh in 1:length(add.range)) {
  for (ii in 1:length(remove.range)) {
    print(paste("clustering. add =", add.range[hh], "remove =", remove.range[ii]))
    
    # get shuffled corum
    tmp.add = addcorum(ints.corum, add.range[hh])
    tmp.remove = removecorum(ints.corum, remove.range[ii])
    # find added + removed interactions
    i.add = !paste(tmp.add$protA, tmp.add$protB, sep="-") %in% paste(ints.corum$protA, ints.corum$protB, sep="-")
    i.remove = !paste(ints.corum$protA, ints.corum$protB, sep="-") %in% paste(tmp.remove$protA, tmp.remove$protB, sep="-")
    # put it together
    ints.noised = rbind(ints.corum[!i.remove,], tmp.add[i.add,])
    
    
    # walktrap
    graph.object = graph_from_edgelist(as.matrix(ints.noised), directed = F)
    walk.cluster = walktrap.community(graph.object)
    for (jj in 1:length(walk.cluster)) {
      if (length(walk.cluster[[jj]]) < 3) next
      cc = cc+1
      clusters$add_mag[cc] = add.range[hh]
      clusters$remove_mag[cc] = remove.range[ii]
      clusters$algorithm[cc] = "walk"
      clusters$cluster[cc] = paste(walk.cluster[[jj]], collapse=";")
    }
    
    # pam
    pam.cluster = pamclust(ints.noised, 100)
    for (jj in 1:length(pam.cluster)) {
      cc = cc+1
      clusters$add_mag[cc] = add.range[hh]
      clusters$remove_mag[cc] = remove.range[ii]
      clusters$algorithm[cc] = "pam"
      clusters$cluster[cc] = pam.cluster[jj]
    }
    
    # mcl
    G = graph.data.frame(ints.noised,directed=FALSE)
    A = as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)
    mcl.cluster = mcl(A, addLoops = FALSE)
    clusts = list()
    unqclusts = unique(mcl.cluster$Cluster)
    for (ii in 1:length(unqclusts)) {
      cc = cc+1
      clusters$add_mag[cc] = add.range[hh]
      clusters$remove_mag[cc] = remove.range[ii]
      clusters$algorithm[cc] = "mcl"
      clusters$cluster[cc] = paste(which(mcl.cluster$Cluster == unqclusts[ii]), collapse = ";")
    }
    
    # clusterone
    co.cluster = clusteroneR(ints.noised, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar")
    for (jj in 1:length(pam.cluster)) {
      clusters$add_mag[cc] = add.range[hh]
      clusters$remove_mag[cc] = remove.range[ii]
      clusters$algorithm[cc] = "co"
      clusters$cluster[cc] = pam.cluster[jj]
    }
  }
}
clusters = clusters[1:cc,]


# calculate Ji1 all clusters (Compare each cluster to its unnoised version)
unqalgs = unique(clusters$algorithm)
clusters$Ji1 = rep(NA, nrow(clusters))
for (kk in 1:length(unqalgs)) {
  I0 = clusters$algorithm==unqalgs[kk]
  ref.clusters = clusters$cluster[I0 & clusters$add_mag==0 & clusters$remove_mag==0]
  for (ii in 1:length(add.range)) {
    ia = clusters$add_mag==add.range[ii]
    for (jj in 1:length(remove.range)) {
      ib = clusters$remove_mag==remove.range[jj]
      I = which(I0 & ia & ib)
      these.clusters = clusters$cluster[I]
      for (mm in 1:length(I)) {
        clusters$Ji1[I[mm]] = calcA(these.clusters[mm], ref.clusters)
      }
    }
  }
}
# write
fn = "../data/clusters_full_netw_walktrap_addANDremove.txt"
write_tsv(clusters, path=fn)
