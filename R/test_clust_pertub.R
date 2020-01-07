
# load clust.perturb
#source("/Users/gregstacey/Academics/Foster/clust-perturb/R/clust-perturb.R")
#source("/Users/gregstacey/Academics/Foster/clust-perturb/R/functions.R")
source("./clust-perturb-tool/clust-perturb.R")
source("./clust-perturb-tool/functions.R")
source("clusterone_java.R")
source("functions.R")

require(igraph)
require(MCL)
require(readr)
require(dplyr)

###################### corum

# binarize human corum
fn = "../data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))
corum = corum[corum$Organism=="Human",]

# turn corum into a list of protein pairs
ints.corum = data.frame(protA = character(10^6),
                        protB = character(10^6), stringsAsFactors = F)
cc = 0
for (ii in 1:nrow(corum)) {
  print(ii)
  prots = sort(unlist(strsplit(corum$`subunits(UniProt IDs)`[ii], ";")))
  if (length(prots)<2) next
  pairs.prots = t(combn(prots, 2))
  
  I = (cc+1) : (cc+nrow(pairs.prots))
  ints.corum$protA[I] = pairs.prots[,1]
  ints.corum$protB[I] = pairs.prots[,2]
  cc = cc+length(I)
}
ints.corum = ints.corum[1:cc,]
ints.corum = distinct(ints.corum)


### cluster
# load("../data/test_clust_pertub.Rda")

noise.range = c(0, 0.01, 0.1, 0.25, 0.5, 0.75)
unqprots = unique(c(ints.corum$protA, ints.corum$protB))
iters = 2
alg.names = c("k-Med", "MCL", "walktrap", "CO")
alg = c(function(x) pam(x, 15),
        function(x) mcl(x, addLoops = FALSE),
        walktrap.community,
        function(x) clusteroneR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar"))
edge.list.format = list(pam.edge.list.format, 
                        mcl.edge.list.format, 
                        function(x) graph_from_edgelist(as.matrix(x), directed = F),
                        NULL)
cluster.format = list(function(x) pam.cluster.format(x,unqprots = unqprots),
                      mcl.cluster.format,
                      NULL,
                      NULL)
clusters.kmed = list()
clusters.mcl = list()
clusters.walk = list()
clusters.co = list()
for (ii in 1:length(noise.range)) {
  print(paste("noise", noise.range[ii]))
  
  # k-med
  print("k-med")
  jj = 1
  clusters.kmed[[ii]] = clust.perturb(ints.corum, clustering.algorithm = alg[[jj]], 
                                      noise = noise.range[ii], iter = iters,
                                      edge.list.format = edge.list.format[[jj]], 
                                      cluster.format = cluster.format[[jj]])
  save(clusters.kmed, clusters.mcl, clusters.walk, clusters.co, 
       file = "../data/test_clust_perturb_4algs.Rda")
  
  # mcl
  print("mcl")
  jj = 2
  clusters.mcl[[ii]] = clust.perturb(ints.corum, clustering.algorithm = alg[[jj]], 
                                     noise = noise.range[ii], iter = iters,
                                     edge.list.format = edge.list.format[[jj]], 
                                     cluster.format = cluster.format[[jj]])
  save(clusters.kmed, clusters.mcl, clusters.walk, clusters.co, 
       file = "../data/test_clust_perturb_4algs.Rda")
  
  # walktrap
  print("walktrap")
  jj = 3
  clusters.walk[[ii]] = clust.perturb(ints.corum, clustering.algorithm = alg[[jj]], 
                                      noise = noise.range[ii], iter = iters,
                                      edge.list.format = edge.list.format[[jj]], 
                                      cluster.format = cluster.format[[jj]])
  save(clusters.kmed, clusters.mcl, clusters.walk, clusters.co, 
       file = "../data/test_clust_perturb_4algs.Rda")
  
  # co
  print("co")
  jj = 4
  clusters.co[[ii]] = clust.perturb(ints.corum, clustering.algorithm = alg[[jj]], 
                                    noise = noise.range[ii], iter = iters,
                                    edge.list.format = edge.list.format[[jj]], 
                                    cluster.format = cluster.format[[jj]])
  save(clusters.kmed, clusters.mcl, clusters.walk, clusters.co, 
       file = "../data/test_clust_perturb_4algs.Rda")
}


### enrichment

# read ontology
ontology = get_ontology("../data/go-basic.obo")
# read annotations
if (1) {
  goa = read_gpa("../data/goa_human.gpa",
                 filter.NOT = T, filter.evidence = c("ND", "IPI", "IEA", "NAS"),
                 ontology = ontology, propagate = T)
  save(goa, file = "../data/goa_human.gpa.Rda")
} else {
  load("../data/goa_human.gpa.Rda")
}
# remove roots 
rootNames = c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
goa %<>% dplyr::filter(!GO.ID %in% rootNames)

# process BP, CC, and MF annotations separately
bp = filter_roots(goa, ontology, 'BP') %>% 
  as_annotation_list("DB_Object_ID", "GO.ID")
cc = filter_roots(goa, ontology, 'CC') %>% 
  as_annotation_list("DB_Object_ID", "GO.ID")
mf = filter_roots(goa, ontology, 'MF') %>% 
  as_annotation_list("DB_Object_ID", "GO.ID")
anns = list(BP = bp, CC = cc, MF = mf)
# create overall annotation object
#ann = as_annotation_list(goa, "DB_Object_ID", "GO.ID")

# filter anns to >5 and <100
unqprots = unique(unlist(ints.corum[,1:2]))
for (ii in 1:length(anns)) {
  print(ii)
  anns[[ii]] = lapply(anns[[ii]], FUN = function(x) { x[x %in% unqprots]
  })
  anns[[ii]] = anns[[ii]][lapply(anns[[ii]],length)>5 & lapply(anns[[ii]],length)<100]
}

# calculate enrichment for every cluster
enr = list()
for (ii in 1:3) {
  enr[[ii]] = data.frame(
    qenriched.goid = character(nrow(clusters)), # which go terms is the clusters enriched for
    np.enriched = numeric(nrow(clusters)), # how many enriched go terms (p<0.01)
    nq.enriched = numeric(nrow(clusters)) # how many enriched go terms (q<0.1)
    , stringsAsFactors = F)
}
for (ii in 1:nrow(clusters)) {
  print(ii)
  for (jj in 1:length(anns)) {
    goids = names(anns[[jj]])
    pp = rep(NA, length(anns[[jj]]))
    for (kk in 1:length(anns[[jj]])) {
      M = sum(unqprots %in% anns[[jj]][[kk]]) # how many proteins have this term (white balls in urn)
      this.cluster = unlist(strsplit(clusters$cluster[ii], ";"))
      K = length(this.cluster) # size of sample (number of balls drawn)
      X = sum(this.cluster %in% anns[[jj]][[kk]]) # go id in this cluster (number of white drawn)
      
      # probability of drawing that many white balls
      pp[kk] = phyper(q=X-1, m=M, n=length(unqprots)-M, k=K, lower.tail=FALSE)
    }
    enr[[jj]]$qenriched.goid[ii] = paste(goids[p.adjust(pp)<.1], collapse = "-")
    enr[[jj]]$np.enriched[ii] = sum(pp<.01)
    enr[[jj]]$nq.enriched[ii] = sum(p.adjust(pp)<.1)
  }
}

summary(glm(enr[[1]]$nq.enriched ~ clusters$reproducibility.J + clusters$size))
summary(glm(enr[[2]]$nq.enriched ~ clusters$reproducibility.J + clusters$size))
summary(glm(enr[[3]]$nq.enriched ~ clusters$reproducibility.J + clusters$size))




###################### our interactomes
# enrichment


