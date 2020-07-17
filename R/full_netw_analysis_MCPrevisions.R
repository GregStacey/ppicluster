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



# command line arguments
this.args = commandArgs(trailingOnly = T)
print(this.args)
param.file = as.character(this.args[1])
param.file.index = as.numeric(this.args[2])

# which system are we on?
home.dir = c("~/projects/ppicluster/",                                     # sockeye
             "/Users/gregstacey/Academics/Foster/ClusterExplore/",         # laptop
             "/home/rstacey/projects/rrg-ljfoster-ab/rstacey/ppicluster/") # cedar
home.dir = home.dir[dir.exists(home.dir)]
R.dir = paste(home.dir, "R/", sep="")
data.dir = paste(home.dir, "data/", sep="")
source(paste(R.dir, "functions.R", sep=""))


# read params
params = as.data.frame(read_tsv(params.file))
# add data.dir to dataset
params$dataset = paste(data.dir, "interactomes/", params$dataset, sep="")

# index params
params = params[params.file.index, ]




# read data
nclust = 1500
data = as.data.frame(read_tsv(params$dataset))
if (grepl("BIOGRID", params$dataset)) {
  data = data[data$`Organism Interactor A`==9606 & data$`Organism Interactor B`==9606, 
              c("Official Symbol Interactor A", "Official Symbol Interactor B")]
} else if (grepl("PE_and_conf_scores", params$dataset)) {
  nclust = 1000
  data = data[order(data$PE_Score, decreasing = T), ][1:9074, c(1, 2, 4)]
} else if (grepl("BioPlex", params$dataset)) {
  data = data[,c("UniprotA", "UniprotB")]
} else if (grepl("email", params$dataset)) {
  nclust = 500
  data = data %>% separate("0 1", c("A", "B")) %>% mutate_all(list(as.numeric=as.numeric)) %>%
    select(c(3, 4)) %>% rename_all(funs(gsub("_as.numeric", "", .))) %>%
    filter(., !A==B)
  for (ii in 1:nrow(data)) data[ii,] = sort(data[ii,])
  data = data[!duplicated(data), ]
  
  data[data==0] = max(data) + 1
}

this.noise.range = params$noise.range
if (this.noise.range == "all") {
  this.noise.range = noise.range
}

for (uu in 1:length(this.noise.range)) {
  sf = paste(data.dir, "clusters/", 
             basename(tools::file_path_sans_ext(params$dataset)),
             "-algorithm=", params$algorithm, 
             "-noise=", this.noise.range[uu], ".txt", sep="")
  print(sf)
  
  # check if output file exists
  if (file.exists(sf)) {
    print(paste("file already exists:", sf))
    #next
  }
  
  # shuffle network
  ints.shuffle = shufflecorum(data, this.noise.range[uu])
  unqprots = unique(c(ints.shuffle[,1], ints.shuffle[,2]))
  if (ncol(ints.shuffle)==2) {names(ints.shuffle) = c("protA", "protB")
  } else names(ints.shuffle) = c("protA", "protB", "weights")
  
  
  # cluster
  if (params$algorithm == "hierarchical") {
    # 1. hierarchical
    x = hierarch.edge.list.format(ints.shuffle)
    tmp = stats::cutree(stats::hclust(d = x, method="average"), k = nclust)
    
    clusts = list()
    for (ii in 1:nclust) {
      clusts[[ii]] = unqprots[tmp == ii]
    }
    
  } else if (params$algorithm == "mcode") {
    # 3. MCODE
    x = graph.data.frame(ints.shuffle)
    #tmp = mcode(x, vwp = 1, haircut = TRUE, fluff = FALSE, fdt = 0.1)
    tmp = mcode(x, vwp = 0, haircut = FALSE, fluff = FALSE, fdt = 0.1)
    clusts = tmp[[1]] %>% lapply(., FUN = function(x) unqprots[x])
    
  } else if (params$algorithm == "louvain") {
    # 4. Louvain
    if (ncol(ints.shuffle)==2) ints.shuffle$weights = 1
    tmp = cluster_resolution(ints.shuffle, 1)
    
    clusts = list()
    unqclusts = unique(tmp$community)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp$community == unqclusts[ii]]
    }
    
  } else if (params$algorithm == "leiden") {
    # 5. Leiden
    adjmat = graph_from_edgelist(as.matrix(ints.shuffle[,1:2]))
    if (ncol(ints.shuffle)==3) edge.attributes(adjmat)$weight = ints.shuffle[,3]
    tmp = leiden(adjmat, resolution_parameter = 0)
    unqprots = V(adjmat)
    
    clusts = list()
    unqclusts = unique(tmp)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp == unqclusts[ii]]
    }
  } else if (params$algorithm == "walk") {
    # walktrap
    graph.object = graph_from_edgelist(as.matrix(ints.shuffle[,1:2]), directed = F)
    if (ncol(ints.shuffle)==3) edge.attributes(graph.object)$weight = ints.shuffle[,3]
    clusts = walktrap.community(graph.object)
    
  } else if (params$algorithm == "pam") {
    # pam (k-med)
    tmp = pamclust(ints.shuffle, nclust)
    clusts = sapply(tmp, strsplit, ";")
    
  } else if (params$algorithm == "co") {
    # cluster one
    tmp = unlist(clusteroneR(ints.shuffle, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar"))
    clusts = sapply(tmp, strsplit, ";")
    
  } else if (params$algorithm == "mcl") {
    # mcl
    G = graph.data.frame(ints.shuffle, directed=FALSE)
    A = as_adjacency_matrix(G,type="both", names=TRUE, sparse=FALSE)
    tmp = MCL::mcl(A, addLoops = FALSE, max.iter = 100)
    unqprots = rownames(A)
    
    clusts = list()
    unqclusts = unique(tmp$Cluster)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp$Cluster == unqclusts[ii]]
    }
  }
  
  # remove clusts with N<3
  nn = unlist(lapply(clusts, length))
  clusts = clusts[nn>2]
  
  
  # write to file
  df = data.frame(cluster = (sapply(clusts, FUN = function(x) paste(x, collapse=";"))),
                  algorithm = rep(params$algorithm, length(clusts)),
                  noise = rep(this.noise.range[uu], length(clusts)),
                  dataset = rep(basename(params$dataset), length(clusts)), stringsAsFactors = F)
  write_tsv(df, path=sf)
}


