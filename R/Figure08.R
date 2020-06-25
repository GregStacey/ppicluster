# predict real world cluster change

source("functions.R")
#source("./clust-perturb-tool/clust-perturb.R")
#source("./clust-perturb-tool/functions.R")
source("clusterone_java.R")

# load  co-fractionation datasets
fn = "../data/interactomes_moredata.txt"
ints.c = as.data.frame(read_tsv(fn))
ints.c = ints.c[ints.c$noise_mag==0, ]
ints.c$ppi = paste(ints.c$protA, ints.c$protB, sep="-")
ints.c = ints.c[ints.c$noise_mag==0,]


## Okay
# after all that, just do replicates!
sets = list(c(1:4),c(5:8), 
            c(9:11),c(12:14), 
            c(15:17),c(18:20), 
            c(21:24),c(25:28))
alg.names = c("k-Med", "MCL", "walktrap", "CO", 
              "hierarchical", "mcode", "louvain", "leiden")
alg = c(function(x) pam(x, 50),
        function(x) mcl(x, infl = 2, remove.self.loops = FALSE),
        walktrap.community,
        function(x) clusteroneR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar"),
        function(x) stats::cutree(stats::hclust(d = hierarch.edge.list.format(x), method="average"), k = 500),
        function(x) mcode(graph.data.frame(x), vwp = 1, haircut = TRUE, fluff = FALSE, fdt = 0.1),
        function(x) {
          x$weights = 1
          return(cluster_resolution(x, 1))},
        function(x) leiden(as_adjacency_matrix(graph_from_edgelist(as.matrix(x))), resolution_parameter = 0))

edge.list.format = list(pam.edge.list.format, 
                        mcl.edge.list.format, 
                        function(x) graph_from_edgelist(as.matrix(x), directed = F),
                        NULL, NULL, NULL, NULL, NULL)
cluster.format = list(pam.cluster.format,
                      mcl.cluster.format,
                      NULL,
                      NULL,
                      hierarch.cluster.format,
                      function(x, unqprots) {    
                        clusts = list()
                        unqclusts = 1:length(x$COMPLEX)
                        for (ii in 1:length(unqclusts)) {
                          clusts[[ii]] = unqprots[x$COMPLEX[[ii]]]
                        }
                        return(clusts)},
                      function(x, unqprots) {    
                        clusts = list()
                        unqclusts = unique(x$community)
                        for (ii in 1:length(unqclusts)) {
                          clusts[[ii]] = unqprots[x$community == unqclusts[ii]]
                        }
                        return(clusts)},
                      function(x, unqprots) {
                        clusts = list()
                        unqclusts = unique(x)
                        for (ii in 1:length(unqclusts)) {
                          tmp = unqprots[x == unqclusts[ii]]
                          if (length(tmp)==0) next
                          clusts[[ii]] = tmp
                        }
                        return(clusts)})
if (F) {
  fns = list.files(path = "../data/", pattern = "^pred_J_expreps_", full.names = T)
  tmp = list()
  for (ii in 1:length(fns)) {
    tmp[[ii]] = as.data.frame(read_tsv(fns[[ii]]))
  }
  df.predrep = bind_rows(tmp, .id = "column_label")
} else {
  # handle command line args
  jj = as.numeric(commandArgs(trailingOnly = T))
  sf = paste("../data/pred_J_expreps_",alg.names[jj],".txt", sep="")
  print(paste("saving data to", sf))
  
  cc = 0
  df.predrep = data.frame(algorithm = character(1e5),
                          predJ = numeric(1e5), # from clust.perturb
                          repJ = numeric(1e5), # comparing to other set
                          set1 = character(1e5),
                          set2 = character(1e5), stringsAsFactors = F)
  for (ii in 1:length(sets)) {
    nn = length(sets[[ii]])
    x1 = sample(sets[[ii]], floor(nn/2))
    x2 = sets[[ii]][!sets[[ii]] %in% x1]
    x1 = sapply(x1, FUN = function(x) paste("us_",x,sep=""))
    x2 = sapply(x2, FUN = function(x) paste("us_",x,sep=""))
    
    # make interactome from subset of replicates
    ia = ints.c$dataset %in% x1 & ints.c$noise_mag==0
    ib = ints.c$dataset %in% x2 & ints.c$noise_mag==0
    protsa = unique(c(ints.c$protA[ia], ints.c$protB[ia]))
    protsb = unique(c(ints.c$protA[ib], ints.c$protB[ib]))
    net1 = ints.c[ia,1:2]
    net2 = ints.c[ib,1:2]
    
    if (nrow(net1)<1000 | nrow(net2)<1000) {
      print(paste("skipping sets", ii))
      next
    }
    
    if (alg.names[jj] == "leiden") {
      x = as.matrix(ints.shuffle)
      adjmat = as_adjacency_matrix(graph_from_edgelist(x))
      tmp = leiden(adjmat, resolution_parameter = 0)
      unqprots = rownames(adjmat)
    } else unqprots = unique(c(net1[,1], net1[,2]))
    
    # cluster with clust.perturb
    #for (jj in 1:4) {
    print(paste("jj=",jj,",  ii=",ii, ", net1"))
    clust1 = clust.perturb(net1, clustering.algorithm = alg[[jj]], noise = 0.15, iters = 2,
                           edge.list.format = edge.list.format[[jj]], 
                           cluster.format = cluster.format[[jj]])
    print(paste("jj=",jj,",  ii=",ii, ", net2"))
    clust2 = clust.perturb(net2, clustering.algorithm = alg[[jj]], noise = 0.15, iters = 2,
                           edge.list.format = edge.list.format[[jj]], 
                           cluster.format = cluster.format[[jj]])
    
    # calc diff b/w sets
    clust1$repJ.clust2 = numeric(nrow(clust1))
    for (ii in 1:nrow(clust1)) {
      clust1$repJ.clust2[ii] = calcA(clust1$cluster[ii], clust2$cluster)
    }
    clust2$repJ.clust1 = numeric(nrow(clust2))
    for (ii in 1:nrow(clust2)) {
      clust2$repJ.clust1[ii] = calcA(clust2$cluster[ii], clust1$cluster)
    }
    
    # store in dataframe
    I = (cc+1) : (cc+nrow(clust1))
    df.predrep$alg[I] = alg.names[jj]
    df.predrep$predJ[I] = clust1$repJ
    df.predrep$repJ[I] = clust1$repJ.clust2
    df.predrep$set1[I] = paste(x1, collapse = ";")
    df.predrep$set2[I] = paste(x2, collapse = ";")
    cc = cc+nrow(clust1)
    I = (cc+1) : (cc+nrow(clust2))
    df.predrep$alg[I] = alg.names[jj]
    df.predrep$predJ[I] = clust2$repJ
    df.predrep$repJ[I] = clust2$repJ.clust1
    df.predrep$set1[I] = paste(x2, collapse = ";")
    df.predrep$set2[I] = paste(x1, collapse = ";")
    cc = cc+nrow(clust2)
    
    # save in case of crash
    write_tsv(df.predrep[1:cc,], path = sf)
  }
}


# read and merge
fns = list.files(path = "../data/", pattern = "^pred_J_expreps_", full.names = T)
tmp = list()
for (ii in 1:length(fns)) {
  tmp[[ii]] = as.data.frame(read_tsv(fns[[ii]]))
}
df.predrep = bind_rows(tmp, .id = "column_label")



ggplot(df.predrep, aes(x=predJ, y=repJ)) + geom_point(alpha = .25) +
  geom_smooth(method = "lm") + facet_grid(~alg) + 
  xlab("Predicted reproducibility (clust.perturb, repJ)") +
  ylab("Actual reproducibility\n(experiment-to-experiment, Ji)") +
  theme_bw()
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/Figure08_v01.pdf"
ggsave(fn, width=10, height=2.6)


