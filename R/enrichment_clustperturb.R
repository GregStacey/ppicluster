source("functions.R")
source("clust-perturb-tool/clust-perturb.R")
source("clust-perturb-tool/functions.R")


# binarize corum
fn = "../data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))
corum = corum[corum$Organism=="Human",]
ints.corum = binarize.corum(corum)
unqprots = unique(c(ints.corum$protA, ints.corum$protB))


#### cluster 4 algorithms

noise = 0.1
iters = 100
alg.names = c("k-Med", "MCL", "walktrap", "CO")
alg = c(function(x) pam(x, 1500),
        function(x) mymcl(x, infl = 2, iter = 100, verbose = T),
        walktrap.community,
        function(x) clusteroneR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar"))
edge.list.format = list(pam.edge.list.format, 
                        mcl.edge.list.format, 
                        function(x) graph_from_edgelist(as.matrix(x), directed = F),
                        NULL)
cluster.format = list(pam.cluster.format,
                      mcl.cluster.format,
                      NULL,
                      NULL)


if (0) {
  load("../data/enrichment_clustperturb.Rda")
} else {
  # # co
  # print("co")
  # jj = 4
  # clusters.co = clust.perturb2(ints.corum, clustering.algorithm = alg[[jj]],
  #                              noise = noise, iter = iters,
  #                              edge.list.format = edge.list.format[[jj]],
  #                              cluster.format = cluster.format[[jj]])
  # save(clusters.kmed, clusters.mcl, clusters.walk, clusters.co,
  #      file = "../data/enrichment_clustperturb.Rda")
  # 
  # # k-med
  # print("k-med")
  # jj = 1
  # clusters.kmed = clust.perturb2(ints.corum, clustering.algorithm = alg[[jj]], 
  #                                noise = noise, iter = iters,
  #                                edge.list.format = edge.list.format[[jj]], 
  #                                cluster.format = cluster.format[[jj]])
  # save(clusters.kmed, clusters.mcl, clusters.walk, clusters.co, 
  #      file = "../data/enrichment_clustperturb.Rda")
  # 
  # # walktrap
  # print("walktrap")
  # jj = 3
  # clusters.walk = clust.perturb2(ints.corum, clustering.algorithm = alg[[jj]],
  #                                noise = noise, iter = iters,
  #                                edge.list.format = edge.list.format[[jj]],
  #                                cluster.format = cluster.format[[jj]])
  # save(clusters.kmed, clusters.mcl, clusters.walk, clusters.co,
  #      file = "../data/enrichment_clustperturb.Rda")

  # mcl
  print("mcl")
  jj = 2
  clusters.mcl = clust.perturb2(ints.corum, clustering.algorithm = alg[[jj]],
                                noise = noise, iter = iters,
                                edge.list.format = edge.list.format[[jj]],
                                cluster.format = cluster.format[[jj]])
  save(clusters.mcl,
       file = "../data/enrichment_clustperturb_mcl.Rda")
}




#### enrichment

# read ontology
ontology = get_ontology("../data/go-basic.obo")
# read annotations
if (0) {
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

# filter anns to >5 and <100
unqprots = unique(unlist(ints.corum[,1:2]))
for (ii in 1:length(anns)) {
  print(ii)
  anns[[ii]] = lapply(anns[[ii]], FUN = function(x) { x[x %in% unqprots]
  })
  anns[[ii]] = anns[[ii]][lapply(anns[[ii]],length)>5 & lapply(anns[[ii]],length)<100]
}


# calculate enrichment for every cluster
clusters = list("clusters.kmed" = clusters.kmed,
                "clusters.mcl" = clusters.mcl,
                "clusters.walk" = clusters.walk,
                "clusters.co" = clusters.co)
all.enr = list() # 4 elements, one for each algorithm
for (aa in c(2,4)) { 
  # loop over algorithms
  
  all.enr[[aa]] = list()
  for (ii in 1:3) {
    all.enr[[aa]][[ii]] = data.frame(
      qenriched.goid = character(nrow(clusters[[aa]])), # which go terms is the clusters[[aa]] enriched for
      np.enriched = numeric(nrow(clusters[[aa]])), # how many enriched go terms (p<0.01)
      nq.enriched = numeric(nrow(clusters[[aa]])) # how many enriched go terms (q<0.1)
      , stringsAsFactors = F)
  }
  for (ii in 1:nrow(clusters[[aa]])) {
    print(ii)
    for (jj in 1:length(anns)) {
      goids = names(anns[[jj]])
      pp = rep(NA, length(anns[[jj]]))
      for (kk in 1:length(anns[[jj]])) {
        M = sum(unqprots %in% anns[[jj]][[kk]]) # how many proteins have this term (white balls in urn)
        this.cluster = unlist(strsplit(clusters[[aa]]$cluster[ii], ";"))
        K = length(this.cluster) # size of sample (number of balls drawn)
        X = sum(this.cluster %in% anns[[jj]][[kk]]) # go id in this cluster (number of white drawn)
        
        # probability of drawing that many white balls
        pp[kk] = phyper(q=X-1, m=M, n=length(unqprots)-M, k=K, lower.tail=FALSE)
      }
      all.enr[[aa]][[jj]]$qenriched.goid[ii] = paste(goids[p.adjust(pp)<.1], collapse = "-")
      all.enr[[aa]][[jj]]$np.enriched[ii] = sum(pp<.01)
      all.enr[[aa]][[jj]]$nq.enriched[ii] = sum(p.adjust(pp)<.1)
    }
  }
  save(clusters, all.enr, file = "../data/enrichment.Rda")
}



