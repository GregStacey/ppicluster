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
this.args = commandArgs(trailingOnly = T)
hparams = -1
if (length(this.args) < 1) {
  print("Testing all sets...")
  hparams = -1
} else if (length(this.args) == 1 ){
  print("Testing one hyperparameter set...")
  hparams = as.integer(as.numeric(this.args))
}



# set up all parameters
# (new data) - (all algorithms)
# (new algorithms) - (all data)
if (dir.exists("/Users/gregstacey/Academics/Foster/Data/dbs/interactomes/")) {
  data.dir = "/Users/gregstacey/Academics/Foster/Data/dbs/interactomes/"
  fns = c("../data/corum_pairwise.txt","../data/ChCh-Miner_durgbank-chem-chem.tsv","../data/email-Eu-core.txt", 
          paste(data.dir,c("BIOGRID-ALL-3.5.186.tab3.txt",
                           "PE_and_conf_scores_01-11-08.txt", # top 9074
                           "BioPlex_293T_Network_10K_Dec_2019.tsv",
                           "HuRI.tsv"), sep=""))

} else {
  data.dir = "../data/interactomes/"
  fns = paste(data.dir,c("corum_pairwise.txt",
                         "ChCh-Miner_durgbank-chem-chem.tsv",
                         "email-Eu-core.txt",
                         "BIOGRID-ALL-3.5.186.tab3.txt",
                         "PE_and_conf_scores_01-11-08.txt", # top 9074
                         "BioPlex_293T_Network_10K_Dec_2019.tsv",
                         "HuRI.tsv"), sep="")
}

algorithms = c("co","pam","mcl","walk", "hierarchical", "mcode", "louvain", "leiden")
noise.range = c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.00)
params = do.call(expand.grid, list(dataset = fns, algorithm = algorithms, noise.range = noise.range)) %>%
  mutate_if(is.factor, as.character) %>% mutate(noise.range = as.numeric(noise.range))
# remove params that are already done
params = params[!((grepl("corum", params$dataset) | grepl("email", params$dataset) | grepl("chem", params$dataset)) &
                     (params$algorithm=="co" | grepl("pam", params$algorithm) | 
                        grepl("mcl", params$algorithm) | grepl("walk", params$algorithm))), ]
# choose which parameter set
if (!hparams==-1) {
  params = params[hparams, ]
}

# check if output file exists
sf = paste("../data/clusters/", 
           basename(tools::file_path_sans_ext(params$dataset)),
           "-algorithm=", params$algorithm, 
           "-noise=", params$noise.range, ".txt", sep="")
if (file.exists(sf)) {
  stop(paste("file already exists:", sf))
} else {
  
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
  
  
  # shuffle network
  ints.shuffle = shufflecorum(data, params$noise.range)
  unqprots = unique(c(ints.shuffle[,1], ints.shuffle[,2]))
  if (ncol(ints.shuffle)==2) {names(ints.shuffle) = c("protA", "protB")
  } else names(ints.shuffle) = c("protA", "protB", "score")
  
  
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
    adjmat = graph_from_edgelist(as.matrix(ints.shuffle))
    tmp = leiden(adjmat, resolution_parameter = 0)
    unqprots = rownames(adjmat)
    
    clusts = list()
    unqclusts = unique(tmp)
    for (ii in 1:length(unqclusts)) {
      clusts[[ii]] = unqprots[tmp == unqclusts[ii]]
    }
  } else if (params$algorithm == "walk") {
    # walktrap
    graph.object = graph_from_edgelist(as.matrix(ints.shuffle), directed = F)
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
    A = as_adjacency_matrix(G,type="both", names=TRUE,sparse=FALSE)
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
                  noise = rep(params$noise.range, length(clusts)),
                  dataset = rep(basename(params$dataset), length(clusts)), stringsAsFactors = F)
  write_tsv(df, path=sf)
}



