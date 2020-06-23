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
    print(paste("walktrap clustering dataset",unqdatasets[jj],"at noise=", unqmags[ii]))
    
    # get shuffled network
    I = ints.c$noise_mag==unqmags[ii] & ints.c$dataset==unqdatasets[jj] & !ints.c$protA==ints.c$protB
    if (sum(I)<190) next
    these.ints = ints.c[I,c("protA", "protB")]
    ints.shuffle = shufflecorum(these.ints, unqmags[ii])

    
    # 3. MCODE
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
    x = ints.shuffle
    x$weights = 1
    tmp = cluster_resolution(x, 1)
    clusts = list()
    unqclusts = unique(tmp$community)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp$community == unqclusts[ii]]
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
    x = as.matrix(ints.shuffle)
    adjmat = as_adjacency_matrix(graph_from_edgelist(x))
    tmp = leiden(adjmat, resolution_parameter = 1)
    clusts = list()
    unqclusts = unique(tmp)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp == unqclusts[ii]]
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
    write_tsv(data.c.acc[1:cc, ], path = "../data/data.c.add_cofracmcp.txt")
  }
}

# write
data.c.add = data.c.add[1:cc,]
write_tsv(data.c.acc[1:cc, ], path = "../data/data.c.add_cofracmcp.txt")

