# Follow reviewer's suggestions to add extra networks and algorithms
# Add these combinations of networksxalgorithms
# (CORUM, emailEU, drugBank) x (hierarchical, DBSCAN, MCODE, Louvain, Leiden)
# (BioGrid, Collins-2007, yeast-consensus, BioPlex, Y2H) x (CO, MCL, k-med, walk, hierarchical, DBSCAN, MCODE, Louvain, Leiden)

# BioGrid = https://downloads.thebiogrid.org/BioGRID
# Collins-2007 = https://pubmed.ncbi.nlm.nih.gov/17200106/
# BioPlex 3.0 = https://www.biorxiv.org/content/10.1101/2020.01.19.905109v1.full.pdf
# Y2H = https://www.nature.com/articles/s41586-020-2188-x (HuRI, http://www.interactome-atlas.org/)
#
# yeast-consensus = https://pubmed.ncbi.nlm.nih.gov/20620961/ (this is just complexes... right?)



source("functions.R")

# command line arguments
hparams = as.integer(as.numeric(commandArgs(trailingOnly = T)))


# set up all parameters
# (new data) - (all algorithms)
# (new algorithms) - (all data)
fns = "../data/interactomes/corum_pairwise.txt"
#algorithms = c("hierarchical", "mcode", "louvain", "leiden")
algorithms = c("hierarchical", "mcode", "louvain", "leiden")
add.range = remove.range = c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.00)
params = do.call(expand.grid, list(dataset = fns, algorithm = algorithms, add_mag = add.range, remove_mag = remove.range)) %>%
  mutate_if(is.factor, as.character) %>% 
  mutate(add_mag = as.numeric(add_mag)) %>% 
  mutate(remove_mag = as.numeric(remove_mag)) 

# choose which parameter set
if (!hparams==-1) {
  params = params[hparams, ]
}

# check if output file exists
sf = paste("../data/clusters/add_remove/", 
           "algorithm=", params$algorithm, 
           "-add=", params$add_mag,
           "-remove=", params$remove_mag,
           ".txt", sep="")
if (file.exists(sf)) {
  stop(paste("file already exists:", sf))
} else {
  
  # read data
  nclust = 1500
  ints.corum = as.data.frame(read_tsv(params$dataset))
  
  
  # shuffle network
  # get shuffled corum
  tmp.add = addcorum(ints.corum, params$add_mag)
  tmp.remove = removecorum(ints.corum, params$remove_mag)
  # find added + removed interactions
  i.add = !paste(tmp.add$protA, tmp.add$protB, sep="-") %in% paste(ints.corum$protA, ints.corum$protB, sep="-")
  i.remove = !paste(ints.corum$protA, ints.corum$protB, sep="-") %in% paste(tmp.remove$protA, tmp.remove$protB, sep="-")
  # put it together
  ints.shuffle = rbind(ints.corum[!i.remove,], tmp.add[i.add,])
  if (nrow(ints.shuffle) < 2500) next  
  
  unqprots = unique(c(ints.shuffle[,1], ints.shuffle[,2]))

  
  # cluster
  if (params$algorithm == "hierarchical") {
    # 1. hierarchical
    unqprots = unique(c(ints.corum$protA, ints.corum$protB))
    x = hierarch.edge.list.format(ints.corum)
    y = stats::hclust(d = x, method="average")
    tmp = stats::cutree(y, k = N.range[ii])
    
    clusts = list()
    for (jj in 1:nclust) {
      clusts[[jj]] = unqprots[tmp == jj]
    }
    
  } else if (params$algorithm == "mcode") {
    # 3. MCODE
    x = graph.data.frame(ints.shuffle)
    clusts = mcode(x, vwp = 1, haircut = TRUE, fluff = FALSE, fdt = 0.1)
    clusts = clusts[[1]] %>% lapply(., FUN = function(x) unqprots[x])
    
  } else if (params$algorithm == "louvain") {
    # 4. Louvain
    x = ints.shuffle
    x$weights = 1
    tmp = cluster_resolution(x, 1)
    
    clusts = list()
    unqclusts = unique(tmp$community)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp$community == unqclusts[ii]]
    }
    
  } else if (params$algorithm == "leiden") {
    # 5. Leiden
    x = as.matrix(ints.shuffle)
    adjmat = as_adjacency_matrix(graph_from_edgelist(x))
    tmp = leiden(adjmat, resolution_parameter = 0)
    unqprots = rownames(adjmat)
    
    clusts = list()
    unqclusts = unique(tmp)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp == unqclusts[ii]]
    }
  } 
  
  # remove clusts with N<3
  nn = unlist(lapply(clusts, length))
  clusts = clusts[nn>2]
  
  
  # write to file
  df = data.frame(cluster = sapply(clusts, FUN = function(x) paste(x, collapse=";")),
                  algorithm = rep(params$algorithm, length(clusts)),
                  add_mag = rep(params$add_mag, length(clusts)),
                  remove_mag = rep(params$remove_mag, length(clusts)), stringsAsFactors = F)
  write_tsv(df, path=sf)
}



