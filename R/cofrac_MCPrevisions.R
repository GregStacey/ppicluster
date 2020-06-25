source("functions.R")

# load clusters
fn1 = "../data/clusters_wshuffle_moredata.txt"
fn2 = "../data/cluster3.txt"
data.c1 = as.data.frame(read_tsv(fn1))
data.c2 = as.data.frame(read_tsv(fn2))
data.c = rbind(data.c1, data.c2)
data.c$noise_mag = as.numeric(data.c$noise_mag)
# reduce to just experiment data
data.c = data.c[!data.c$data_type%in%"corum",]

# load interactomes
fn = "../data/interactomes_moredata.txt"
ints.c = as.data.frame(read_tsv(fn))
ints.c$ppi = paste(ints.c$protA, ints.c$protB, sep="-")

unqalgs = unique(data.c$algorithm)
unqdatasets = unique(data.c$data_type)
unqmags = sort(unique(data.c$noise_mag))
unqmags = unqmags[unqmags<=1]

# walktrap cluster
# add 1 iterations of mcode, louvain, and leiden
data.c.add = data.frame(data_type = character(10^5), noise_type = character(10^5),
                        noise_mag = numeric(10^5), algorithm = character(10^5),
                        cluster = character(10^5), stringsAsFactors = F)
cc = 0
for (ii in 1:length(unqmags)) {
  for (jj in 1:length(unqdatasets)) {
    print(paste("clustering dataset",unqdatasets[jj],"at noise=", unqmags[ii]))
    
    # get shuffled network
    I = ints.c$noise_mag==unqmags[ii] & ints.c$dataset==unqdatasets[jj] & !ints.c$protA==ints.c$protB
    if (sum(I, na.rm=T)<190) next
    these.ints = ints.c[I,c("protA", "protB")]
    ints.shuffle = shufflecorum(these.ints, unqmags[ii])
    unqprots = unique(c(ints.shuffle[,1], ints.shuffle[,2]))
    
    
    # 3. MCODE
    print("mcode")
    x = graph.data.frame(ints.shuffle)
    clusts = mcode(x, vwp = 1, haircut = TRUE, fluff = FALSE, fdt = 0.1)
    clusts = clusts[[1]] %>% lapply(., FUN = function(x) unqprots[x])
    for (kk in 1:length(clusts)) {
      if (length(clusts[[kk]]) < 3) next
      cc = cc+1
      data.c.add$data_type[cc] = unqdatasets[jj]
      data.c.add$noise_type[cc] = "chrom"
      data.c.add$noise_mag[cc] = unqmags[ii]
      data.c.add$algorithm[cc] = "mcode"
      data.c.add$cluster[cc] = paste(clusts[[kk]], collapse=";")
    }
    
    
    # 4. Louvain
    print("louvain")
    x = ints.shuffle
    x$weights = 1
    tmp = cluster_resolution(x, 1)
    clusts = list()
    unqclusts = unique(tmp$community)
    for (uu in 1:length(unqclusts)) {
      clusts[[uu]] = unqprots[tmp$community == unqclusts[uu]]
    }
    for (kk in 1:length(clusts)) {
      if (length(clusts[[kk]]) < 3) next
      cc = cc+1
      data.c.add$data_type[cc] = unqdatasets[jj]
      data.c.add$noise_type[cc] = "chrom"
      data.c.add$noise_mag[cc] = unqmags[ii]
      data.c.add$algorithm[cc] = "louvain"
      data.c.add$cluster[cc] = paste(clusts[[kk]], collapse=";")
    }
    
    
    # 5. Leiden
    print("leiden")
    # 5. Leiden
    x = as.matrix(ints.shuffle)
    adjmat = as_adjacency_matrix(graph_from_edgelist(x))
    tmp = leiden(adjmat, resolution_parameter = 1)
    unqprots = rownames(adjmat)
    clusts = list()
    unqclusts = unique(tmp)
    for (jj in 1:length(unqclusts)) {
      clusts[[jj]] = unqprots[tmp == unqclusts[jj]]
    }
    for (kk in 1:length(clusts)) {
      if (length(clusts[[kk]]) < 3) next
      cc = cc+1
      data.c.add$data_type[cc] = unqdatasets[jj]
      data.c.add$noise_type[cc] = "chrom"
      data.c.add$noise_mag[cc] = unqmags[ii]
      data.c.add$algorithm[cc] = "leiden"
      data.c.add$cluster[cc] = paste(clusts[[kk]], collapse=";")
    }
    
    # write in case of crash
    write_tsv(data.c.add[1:cc, ], path = "../data/data.c.add_cofracmcp.txt")
  }
}

# write
data.c.add = data.c.add[1:cc,]
write_tsv(data.c.add[1:cc, ], path = "../data/data.c.add_cofracmcp.txt")



# calculate Ji
nn = 10^6
df = data.frame(dataset = character(nn), algorithm = character(nn), noise_mag=numeric(nn),
                nclust0 = numeric(nn), nclust = numeric(nn), clust.size=numeric(nn), clustJ=numeric(nn),
                nint0 = numeric(nn), nint = numeric(nn), intJ=numeric(nn),
                stringsAsFactors = F)
unqdatasets = unique(data.c.add$data_type)
unqmags = unique(data.c.add$noise_mag)
unqalgs = unique(data.c.add$algorithm)
cc = 0
for (ii in 1:length(unqdatasets)) {
  print(paste("   ", ii))
  i0 = ints.c$dataset%in%unqdatasets[ii] & ints.c$noise_mag==0
  this.nint0 = sum(i0)
  for (jj in 1:length(unqmags)) {
    print(unqmags[jj])
    i1 = ints.c$dataset%in%unqdatasets[ii] & ints.c$noise_mag==unqmags[jj]
    this.intJ = length(intersect(ints.c$ppi[i0], ints.c$ppi[i1])) / 
      length(unique(c(ints.c$ppi[i0],ints.c$ppi[i1])))
    this.nint = sum(i1)
    for (kk in 1:length(unqalgs)) {
      I0 = data.c.add$data_type%in%unqdatasets[ii] & data.c.add$noise_mag==0 & data.c.add$algorithm%in%unqalgs[kk]
      I1 = data.c.add$data_type%in%unqdatasets[ii] & data.c.add$noise_mag==unqmags[jj] & data.c.add$algorithm%in%unqalgs[kk]
      if (sum(I1)==0) next
      
      cluster0 = data.c.add$cluster[I0]
      cluster = data.c.add$cluster[I1]
      J.i = numeric(length(cluster))
      for (mm in 1:length(cluster)) {
        J.i[mm] = calcA(cluster[mm], cluster0)
        
        cc = cc+1
        df$dataset[cc] = unqdatasets[ii]
        df$algorithm[cc] = unqalgs[kk]
        df$noise_mag[cc] = unqmags[jj]
        df$nclust0[cc] = length(cluster0)
        df$nclust[cc] = length(cluster)
        df$clustJ[cc] = J.i[mm]
        df$nint0[cc] = this.nint0
        df$nint[cc] = this.nint
        df$intJ[cc] = this.intJ
        df$clust.size[cc] = length(unlist(strsplit(cluster[mm], ";")))
      }
    }
  }
}
df = df[1:cc,]

# write finally
write_tsv(cbind(data.c.add, df), path = "../data/data.c.add_cofracmcp.txt")
