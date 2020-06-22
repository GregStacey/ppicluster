# calculate Ji for all algorithm/datasets

if (dir.exists("/Users/gregstacey/Academics/Foster/Data/dbs/interactomes/")) {
  data.dir = "/Users/gregstacey/Academics/Foster/Data/dbs/interactomes/"
  fns = paste(data.dir,c("BIOGRID-ALL-3.5.186.tab3.txt",
                         "PE_and_conf_scores_01-11-08.txt", # top 9074
                         "BioPlex_293T_Network_10K_Dec_2019.tsv",
                         "HuRI.tsv"), sep="")
  fns = c("../data/corum_pairwise.txt","../data/ChCh-Miner_durgbank-chem-chem.tsv","../data/email-Eu-core.txt", 
          fns)
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
datasets = unique(basename(tools::file_path_sans_ext(fns)))
params2 = do.call(expand.grid, list(dataset = fns, algorithm = algorithms))

noise.range = c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.00)

# check what files exist
ii = as.integer(as.numeric(commandArgs(trailingOnly = T)))
# what noise ranges are there?
sfs = paste("../data/clusters/", 
            basename(tools::file_path_sans_ext(params2$dataset[ii])),
            "-algorithm=", params2$algorithm[ii], 
            "-noise=", noise.range, ".txt", sep="")

# if there's not a full noise set, abort
if (!sum(sapply(sfs, file.exists)) == length(noise.range)){
  stop(paste("full noise range not run for dataset", 
             params2$dataset[ii],
             "and algorithm", params2$algorithm[ii]))
}


# read clusters
clusts = as.data.frame(read_tsv(sfs[1])) 
names(clusts) = c("cluster", "algorithm", "noise_mag", "network")
clusters = clusts[,c(4,3,2,1)]
clusters$Ji = 1
for (jj in 2:length(noise.range)) {
  # read and add clusters
  clusts = as.data.frame(read_tsv(sfs[jj]))[, c(4,3,2,1)] 
  names(clusts) = c("network", "noise_mag", "algorithm", "cluster")
  clusts$Ji = rep(NA, nrow(clusts))
  clusters = rbind(clusters, clusts[,])

  # calculate Ji
  I0 = which(clusters$noise_mag==0)
  I = which(clusters$noise_mag==noise.range[jj])
  these.clusters = clusters$cluster[I]
  ref.clusters = clusters$cluster[I0]
  for (mm in 1:length(I)) {
    clusters$Ji[I[mm]] = calcA(these.clusters[mm], ref.clusters)
  }
}

sf = paste("../data/clusters/outcomes/", 
           basename(tools::file_path_sans_ext(params2$dataset[ii])),
           "-algorithm=", params2$algorithm[ii], ".txt", sep="")
write_tsv(clusters, path=sf)

