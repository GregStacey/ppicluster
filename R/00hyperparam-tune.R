# maximize silhouette score for CORUM and these algorithms
# - hierarchical
# - DBSCAN
# - MCODE
# - Louvain
# - Leiden

  
source("functions.R")


# command line arguments
this.args = commandArgs(trailingOnly = T)
hparams = -1
if (length(this.args) < 1) {
  print("Testing all sets...")
  hparams = -1
} else if (length(this.args) == 1 ){
  print("Testing one hyperparameter set...")
  hparams = as.integer(as.numeric(this.args))
}


# binarize corum
fn = "../data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))
corum = corum[corum$Organism=="Human",]

# turn corum into a list of protein pairs
ints.corum = as.data.frame(read_tsv("../data/interactomes/corum_pairwise.txt"))

unqprots = unique(c(ints.corum[,1], ints.corum[,2]))



# set up all paramaters
df.params = data.frame(algorithm = character(1e4),
                       params = character(1e4), stringsAsFactors = F)
cc = 0
# hierarchical
nclusters = c(50, 100, 250, 500, 1000, 1500, 2000, 5000)
for (ii in 1:length(nclusters)) {
  cc = cc+1
  df.params[cc,] = c("hierarchical", nclusters[ii])
}
# dbscan
eps.range = c(.01, .1, .2, .5, 1, 2)
minPoints = c(2, 3, 5, 10)
tmp = unlist(apply(do.call(expand.grid, list(eps.range = eps.range, minPoints = minPoints)), 
                   1, paste, collapse = " - "))
for (ii in 1:length(tmp)) {
  cc = cc+1
  df.params[cc,] = c("dbscan", tmp[ii])
}
# mcode
haircut = c(T, F)
fluff = c(T, F)
fluff.range = c(0.1, 0.3, 0.5, 0.7, 0.9)
vwp.range = c(0.1, 0.3, 0.5, 0.7, 0.9)
tmp = unlist(apply(do.call(expand.grid, list(haircut = haircut, fluff = fluff, fluff.range = fluff.range, vwp.range = vwp.range)), 
                   1, paste, collapse = " - "))
for (ii in 1:length(tmp)) {
  cc = cc+1
  df.params[cc,] = c("mcode", tmp[ii])
}
# louvain
res.range = c(0, .25, 0.5, .75, 1, 1.5, 2, 4, 8, 10, 15, 20)
for (ii in 1:length(res.range)) {
  cc = cc+1
  df.params[cc,] = c("louvain", res.range[ii])
}
# leiden
res.range = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
for (ii in 1:length(res.range)) {
  cc = cc+1
  df.params[cc,] = c("leiden", res.range[ii])
}
df.params = df.params[1:cc, ]

# hack hack hack hack hack hack
# reduce to just louvain
df.params = df.params[df.params$algorithm=="louvain", ]
#

df.params0 = df.params
for (hparams in 1:nrow(df.params0)) {
  # choose which parameter set
  if (!hparams==-1) {
    df.params = df.params0[hparams, ]
    params = as.numeric(unlist(strsplit(df.params$params, " - ")))
  }

  # cluster
  if (df.params$algorithm == "hierarchical") {
    # 1. hierarchical
    x = hierarch.edge.list.format(ints.corum)
    tmp = stats::cutree(stats::hclust(d = x, method="average"), k = params[1])
    
    clusts = list()
    for (ii in 1:params[1]) {
      clusts[[ii]] = unqprots[tmp == ii]
    }
    
  } else if (df.params$algorithm == "dbscan") {
    # 2. DBSCAN
    # https://towardsdatascience.com/how-dbscan-works-and-why-should-i-use-it-443b4a191c80
    eps.range = c(.01, .1, .2, .5, 1, 2)
    minPoints = c(2, 3, 5, 10)
    x = get.adjacency(graph.data.frame(ints.corum))
    tmp = dbscan(x, eps = params[1], params[2])
    
    clusts = list()
    unqclusts = unique(tmp$cluster)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp$cluster == unqclusts[ii]]
    }
    
  } else if (df.params$algorithm == "mcode") {
    # 3. MCODE
    x = graph.data.frame(ints.corum)
    clusts = mcode(x, vwp = params[4], haircut = as.logical(params[1]), fluff = as.logical(params[2]), fdt = params[3])
    clusts = clusts[[1]] %>% lapply(., FUN = function(x) unqprots[x])
    
  } else if (df.params$algorithm == "louvain") {
    # # 4. Louvain
    # tmp = louvain(ints.corum, params[1])
    # 
    # clusts = list()
    # unqclusts = unique(tmp$cluster)
    # for (ii in 1:length(unqclusts)) {
    #   clusts[[ii]] = unqprots[tmp$clust == unqclusts[ii]]
    # }
    
    # 4. Louvain
    if (ncol(ints.corum)==2) ints.corum$weights = 1
    tmp = louvain(ints.corum, params[1])
    
    clusts = list()
    unqclusts = unique(tmp$cluster)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp$cluster == unqclusts[ii]]
    }
    # remove singletons
    clusts = clusts[!tolower(unqclusts) =="singleton"]
    
  } else if (df.params$algorithm == "leiden") {
    # # 5. Leiden
    # adjmat = graph_from_edgelist(as.matrix(ints.corum))
    # tmp = leiden(adjmat, resolution_parameter = params[1])
    # unqprots = rownames(adjmat)
    # 
    # clusts = list()
    # unqclusts = unique(tmp)
    # for (ii in 1:length(unqclusts)) {
    #   clusts[[ii]] = unqprots[tmp == unqclusts[ii]]
    # }
    
    # 5. Leiden
    adjmat = graph_from_edgelist(as.matrix(ints.corum[,1:2]))
    if (ncol(ints.corum)==3) edge.attributes(adjmat)$weight = ints.corum[,3]
    this.clust = leiden(adjmat, resolution_parameter = params[1])
    unqprots = names(V(adjmat))
    
    unqclusts = sort(unique(this.clust))
    clusts = list()
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[this.clust == unqclusts[ii]]
    }
  }
  
  # remove clusts with N<3
  nn = unlist(lapply(clusts, length))
  clusts = clusts[nn>2]

  # compare to ground truth
  df = data.frame(cluster = unlist(lapply(clusts, paste, collapse = ";")), stringsAsFactors = F)
  df$Ji = unlist(sapply(df$cluster, FUN = function(x) calcA(x, corum$`subunits(UniProt IDs)`)))
  df$algorithm = df.params$algorithm
  df$params = df.params$params
  
  # write
  sf = paste("../data/00hyperparam_tune/louvain_", hparams, "_paramtune.txt", sep="")
  write_tsv(df, path = sf)  
}



# read all files and get optimal sets for each algorithm
fns = list.files("../data/00hyperparam_tune/", pattern = "louvain*", full.names = T)
df = as.data.frame(read_tsv(fns[1]))
for (ii in 1:length(fns)) {
  print(ii)
  df = rbind(df, as.data.frame(read_tsv(fns[ii])))
}
df = df[!is.na(df$algorithm), ]

df.best = data.frame(algorithm = unique(df$algorithm),
                     J = rep(NA, length(unique(df$algorithm))),
                     params = rep(NA, length(unique(df$algorithm))), stringsAsFactors = F)
for (ii in 1:nrow(df.best)) {
  I = df$algorithm == df.best$algorithm[ii]
  unqparams = sort(unique(df$params[I]))
  JJ = numeric(length(unqparams))
  for (jj in 1:length(unqparams)) {
    ia = I & df$params == unqparams[jj]
    JJ[jj] = median(df$Ji[ia], na.rm=T)
  }
  
  ib = which.max(JJ)
  df.best$params[ii] = unqparams[ib]
  df.best$J[ii] = max(JJ, na.rm=T)
  print(paste(df.best$algorithm[ii], df.best$params[ii], sep=" : "))
}
write_tsv(df.best, path = "../data/00hyperparam_tune/optim_params.txt")




## debug
# do some algorithms really cluster so poorly?

# 1. hierarchical
x = as.dist(get.adjacency(graph.data.frame(ints.corum)))
tmp = stats::cutree(stats::hclust(d = x, method="complete"), k = 500)
clusts = list()
for (ii in 1:500) {
  clusts[[ii]] = unqprots[tmp == ii]
}

meths = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
nn1 = nn = numeric(length(meths))
for (ii in 1:length(meths)) {
  print(paste(ii, meths[ii]))
  tmp = stats::cutree(stats::hclust(d = 1-x, method=meths[ii]), k = 2000)
  nn[ii] = sd(table(tmp))
  
  x1 = pam.edge.list.format(ints.corum)
  tmp1 = stats::cutree(stats::hclust(d = x1, method=meths[ii]), k = 2000)
  nn1[ii] = sd(table(tmp1))
}

df = data.frame(cluster = unlist(lapply(clusts, paste, collapse = ";")), stringsAsFactors = F)
df$Ji = unlist(sapply(df$cluster, FUN = function(x) calcA(x, corum$`subunits(UniProt IDs)`)))
df$algorithm = df.params$algorithm
df$params = df.params$params
df$nsize = unlist(sapply(df$cluster, FUN = function(x) length(unlist(strsplit(x, ";")))))


