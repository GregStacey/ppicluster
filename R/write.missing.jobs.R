

# which system are we on?
home.dir = c("~/projects/ppicluster/",                                     # sockeye
             "/Users/gregstacey/Academics/Foster/ClusterExplore/",         # laptop
             "/home/rstacey/projects/rrg-ljfoster-ab/rstacey/ppicluster/") # cedar
home.dir = home.dir[dir.exists(home.dir)]
R.dir = paste(home.dir, "R/", sep="")
data.dir = paste(home.dir, "data/", sep="")

source("functions.R")


# which jobs are missing
algorithms = c("co","pam","mcl","walk", "hierarchical", "mcode", "louvain", "leiden")
noise.range = c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.00)
dbs = c("corum_pairwise.txt", "ChCh-Miner_durgbank-chem-chem.tsv", "email-Eu-core.txt",
        "BIOGRID-ALL-3.5.186.tab3.txt", "PE_and_conf_scores_01-11-08.txt", 
        "BioPlex_293T_Network_10K_Dec_2019.tsv", "HuRI.tsv")
jobs = missing.jobs(paste(data.dir, "/clusters", sep=""), algs = algorithms,
                    dbs = dbs, noise.range = noise.range)


# write missing jobs
# ia = ...
write_tsv(jobs[ia,], path = "../data/jobs_biogrid_co.txt")
