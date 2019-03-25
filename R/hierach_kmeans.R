# Add extra clustering algorithms to MCL, CL1, and MCL+CL1

source("functions.R")

# load interactomes to cluster
fn1 = "../data/interactomes_moredata.txt"
fn2 = "../data/interactomes_corum_shuffle.txt"
ints1 = as.data.frame(read_tsv(fn1))
ints2 = as.data.frame(read_tsv(fn2))
ints.c = rbind(ints1, ints2)
ints.c$ppi = paste(ints.c$protA, ints.c$protB, sep="-")
# remove self-interactions (caused by removing isoform tags)
ints.c = ints.c[!ints.c$protA==ints.c$protB,]
unqdatasets = sort(unique(ints.c$dataset))
unqnoise = unique(ints.c$noise_mag)

# load previous clusters
fn = "../data/clusters_wshuffle_moredata.txt"
clusters = as.data.frame(read_tsv(fn))
clusters$noise_mag = as.numeric(clusters$noise_mag)

# make a duplicate dataframe to fill with new clusters
clusters2 = data.frame(data_type = character(10^6),
                       noise_type = character(10^6),
                       noise_mag = numeric(10^6),
                       algorithm = character(10^6),
                       cluster = character(10^6), stringsAsFactors = F)




# hierarchical cluster
cc = 0
for (ii in 1:length(unqdatasets)) {
  nclust = 400
  if (unqdatasets[ii]=="corum") {
    nclust = 1500
  }
  for (jj in 1:length(unqnoise)) {
    print(paste(unqdatasets[ii], unqnoise[jj]))
    
    # cluster
    I = ints.c$dataset==unqdatasets[ii] & ints.c$noise_mag==unqnoise[jj]
    ints = ints.c[I,]
    ints = ints[!duplicated(ints$ppi),]
    if (nrow(ints)<nclust | length(unique(c(ints$protA, ints$protB)))<nclust) {
      print("   skipping")
      next
    }
    clusts = hiclust(ints,nclust)
    
    # quality control
    if (TRUE) { # confirm that clusters have interactions in them!!
      dd = numeric(length(clusts))
      for (kk in 1:length(clusts)) {
        prots = unlist(strsplit(clusts[kk], ";"))
        #prots = sample(unique(c(ints$protA, ints$protA)), length(prots))
        nn = length(prots)
        ia = ints$protA%in%prots & ints$protB%in%prots
        dd[kk] = sum(ia) / (nn * (nn-1) / 2)
      }
      
      if (sum(dd==0)>0) {
        error("Awooga!!")
      }
    }
    
    # append to dataframe
    I.append = (cc+1) : (cc+length(clusts))
    clusters2$data_type[I.append] = unqdatasets[ii]
    clusters2$noise_type[I.append] = "chrom"
    if (unqdatasets[ii]=="corum") clusters2$noise_type[I.append] = "network_shuffle"
    clusters2$noise_mag[I.append] = unqnoise[jj]
    clusters2$algorithm[I.append] = "hierarchical"
    clusters2$cluster[I.append] = clusts
    
    cc = cc+length(I.append)
    print(length(clusts))
  }
}



# PAM cluster
for (ii in 1:length(unqdatasets)) {
  nclust = 400
  if (unqdatasets[ii]=="corum") {
    nclust = 900
  }
  for (jj in 1:length(unqnoise)) {
    print(paste(unqdatasets[ii], unqnoise[jj]))
    
    # cluster
    I = ints.c$dataset==unqdatasets[ii] & ints.c$noise_mag==unqnoise[jj]
    ints = ints.c[I,]
    ints = ints[!duplicated(ints$ppi),]
    if (nrow(ints)<nclust | length(unique(c(ints$protA, ints$protB)))<nclust) {
      print("   skipping")
      next
    }
    clusts = pamclust(ints,nclust)
    
    
    # quality control
    if (TRUE) { # confirm that clusters have interactions in them!!
      dd = numeric(length(clusts))
      for (kk in 1:length(clusts)) {
        prots = unlist(strsplit(clusts[kk], ";"))
        #prots = sample(unique(c(ints$protA, ints$protA)), length(prots))
        nn = length(prots)
        ia = ints$protA%in%prots & ints$protB%in%prots
        dd[kk] = sum(ia) / (nn * (nn-1) / 2)
      }
      
      if (sum(dd==0)>0) {
        error("Awooga!!")
      }
    }
    
    # append to dataframe
    I.append = (cc+1) : (cc+length(clusts))
    clusters2$data_type[I.append] = unqdatasets[ii]
    clusters2$noise_type[I.append] = "chrom"
    if (unqdatasets[ii]=="corum") clusters2$noise_type[I.append] = "network_shuffle"
    clusters2$noise_mag[I.append] = unqnoise[jj]
    clusters2$algorithm[I.append] = "pam"
    clusters2$cluster[I.append] = clusts
    
    cc = cc+length(I.append)
    print(length(clusts))
  }
}




# write clusters2
clusters2 = clusters2[1:cc,]
fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/cluster2.txt"
write_tsv(clusters2, fn)

