
### motivate the tool
# show there are clusters that just tend to break.

source("functions.R")

# get binarized corum
# load("../data/corum_old/all_corum_ints.Rda")
# ints.corum = ints.old.corum[[5]]
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
unqprots = unique(c(ints.corum$protA, ints.corum$protB))

# cluster corum
co.alg = function(x) clusteroneR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar")
clusts = unlist(co.alg(ints.corum))
nn = length(clusts)

# binarize the clusters, assign each edge to a cluster
df = data.frame(protA = character(1e6), protB = character(1e6),
                cluster = numeric(1e6), stringsAsFactors = F)
cc = 0
for (ii in 1:length(clusts)) {
  prots = sort(unlist(strsplit(clusts[ii], ";")))
  if (length(prots)<2) next
  pairs.prots = t(combn(prots, 2))
  
  I = (cc+1) : (cc+nrow(pairs.prots))
  df$protA[I] = pairs.prots[,1]
  df$protB[I] = pairs.prots[,2]
  df$cluster[I] = ii
  cc = cc+length(I)
}
df = df[1:cc, ]

iterMax = 50
nn2 = iterMax * length(clusts)
cc = 0
broken = data.frame(id.clust = numeric(nn2),
                    iter = numeric(nn2),
                    Ji = numeric(nn2),
                    n.scrambled = numeric(nn2),
                    frac.scrambled = numeric(nn2), stringsAsFactors = F)
for (iter in 1:iterMax) {
  print(paste("iter =", iter))
  # pick 1/4 of the clusters
  scrambledclusts = sample(nn, round(nn/4))
  # define nodes as i) scrambled, ii) connected to a scrambled node, or iii) segregated
  I.scrambled = df$cluster %in% scrambledclusts
  scrambledprots = unique(c(df$protA[I.scrambled==1], df$protB[I.scrambled==1]))
  segregatedprots = unqprots[!unqprots %in% scrambledprots]
  
  df$connected = df$protA %in% scrambledprots | df$protB %in% scrambledprots
  df$segregated = !I.scrambled & !df$connected
  segregatedprots = unique(c(df$protA[df$segregated==1], df$protB[df$segregated==1]))
  
  # identify which edges to shuffle in ints.corum
  # i.e. between two scrambled prots
  I = ints.corum$protA %in% scrambledprots & ints.corum$protB %in% scrambledprots
  # shuffle these edges with a fairly high noise level
  shuffled.part = shufflecorum(ints.corum[I,], 0.75)
  # combine
  shuffled.ints = rbind(ints.corum[!I,], shuffled.part)
  
  # confirm all edges between segregated prots are unchanged
  ia = ints.corum$protA %in% segregatedprots & ints.corum$protB %in% segregatedprots
  ints0 = sort(paste(ints.corum$protA[ia], ints.corum$protB[ia], sep="-"))
  ib = shuffled.ints$protA %in% segregatedprots & shuffled.ints$protB %in% segregatedprots
  ints1 = sort(paste(shuffled.ints$protA[ib], shuffled.ints$protB[ib], sep="-"))
  if (!identical(ints0, ints1)) error
  
  # cluster shuffled.ints
  clusts2 = unlist(co.alg(shuffled.ints))
  
  # calculate Ji
  for (jj in 1:length(clusts)) {
    cc = cc+1
    prots = unlist(strsplit(clusts[jj], ";"))
    broken$id.clust[cc] = jj
    broken$iter[cc] = iter
    broken$Ji[cc] = calcA(clusts[jj], clusts2)
    broken$n.scrambled[cc] = sum(prots %in% scrambledprots)
    broken$frac.scrambled[cc] = broken$n.scrambled[cc] / length(prots)
  }
  save(ints.corum, clusts, broken, file="../data/broken.Rda")
}
save(ints.corum, clusts, broken, file="../data/broken.Rda")




# # test each cluster
# #   Ji ~ f(n.scrambled) + g(cluster.id)
# #   is g significant
# #   only use data from clusters with size +/- nn
# unique(broken$id.clust[broken$n.scrambled==0 & broken$Ji<0.5])
# 
# clusters = data.frame(clust = clusts,
#                       size = numeric(length(clusts)),
#                       pp = numeric(length(clusts)),
#                       effect = numeric(length(clusts)), stringsAsFactors = F)
# for (ii in 1:length(clusts)) {
#   I = broken$id.clust == ii
#   nn = unique(broken$size[I])
#   
#   # find at least 10 other clusters with similar size
#   cc = 0
#   I.background = broken$size %in% c((nn-cc):(nn+cc))
#   while (sum(I.background) < (sum(I)*10)) {
#     cc = cc+1
#     I.background = broken$size %in% c((nn-cc):(nn+cc))
#   }
#   
#   df = data.frame(Ji = broken$Ji[I.background],
#                   n.scrambled = broken$n.scrambled[I.background],
#                   frac.scrambled = broken$frac.scrambled[I.background],
#                   id = broken$id.clust[I.background]==ii)
#   #tmp = summary(glm(Ji ~ n.scrambled + id , data=df))
#   tmp = summary(glm(Ji ~ n.scrambled + id + id*n.scrambled, data=df))
#   clusters$pp[ii] = tmp$coefficients[rownames(tmp$coefficients)=="idTRUE", 4]
#   clusters$size[ii] = nn
#   clusters$effect[ii] = tmp$coefficients[rownames(tmp$coefficients)=="idTRUE", 1]
#   
#     #ggplot(df, aes(x=n.scrambled, y=Ji)) + geom_point(alpha=.25) + facet_wrap(~id)
# }
# clusters$qq = p.adjust(clusters$pp)
# ggplot(clusters, aes(x=abs(effect), y=-log10(qq), color=effect>0)) + geom_point(alpha=.5)
# 
# 
# 
# # does pp correlate with enrichment?
# # read ontology
# ontology = get_ontology("../data/go-basic.obo")
# # read annotations
# if (0) {
#   goa = read_gpa("../data/goa_human.gpa",
#                  filter.NOT = T, filter.evidence = c("ND", "IPI", "IEA", "NAS"),
#                  ontology = ontology, propagate = T)
#   save(goa, file = "../data/goa_human.gpa.Rda")
# } else {
#   load("../data/goa_human.gpa.Rda")
# }
# # remove roots 
# rootNames = c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
# goa %<>% dplyr::filter(!GO.ID %in% rootNames)
# 
# # process BP, CC, and MF annotations separately
# bp = filter_roots(goa, ontology, 'BP') %>% 
#   as_annotation_list("DB_Object_ID", "GO.ID")
# cc = filter_roots(goa, ontology, 'CC') %>% 
#   as_annotation_list("DB_Object_ID", "GO.ID")
# mf = filter_roots(goa, ontology, 'MF') %>% 
#   as_annotation_list("DB_Object_ID", "GO.ID")
# anns = list(BP = bp, CC = cc, MF = mf)
# # create overall annotation object
# #ann = as_annotation_list(goa, "DB_Object_ID", "GO.ID")
# 
# # filter anns to >5 and <100
# unqprots = unique(unlist(ints.corum[,1:2]))
# for (ii in 1:length(anns)) {
#   print(ii)
#   anns[[ii]] = lapply(anns[[ii]], FUN = function(x) { x[x %in% unqprots]
#   })
#   anns[[ii]] = anns[[ii]][lapply(anns[[ii]],length)>5 & lapply(anns[[ii]],length)<100]
# }
# 
# # calculate enrichment for every cluster
# enr = list()
# for (ii in 1:3) {
#   enr[[ii]] = data.frame(
#     qenriched.goid = character(nrow(clusters)), # which go terms is the clusters enriched for
#     np.enriched = numeric(nrow(clusters)), # how many enriched go terms (p<0.01)
#     nq.enriched = numeric(nrow(clusters)) # how many enriched go terms (q<0.1)
#     , stringsAsFactors = F)
# }
# for (ii in 1:nrow(clusters)) {
#   print(ii)
#   for (jj in 1:length(anns)) {
#     goids = names(anns[[jj]])
#     pp = rep(NA, length(anns[[jj]]))
#     for (kk in 1:length(anns[[jj]])) {
#       M = sum(unqprots %in% anns[[jj]][[kk]]) # how many proteins have this term (white balls in urn)
#       this.cluster = unlist(strsplit(clusters$clust[ii], ";"))
#       K = length(this.cluster) # size of sample (number of balls drawn)
#       X = sum(this.cluster %in% anns[[jj]][[kk]]) # go id in this cluster (number of white drawn)
#       
#       # probability of drawing that many white balls
#       pp[kk] = phyper(q=X-1, m=M, n=length(unqprots)-M, k=K, lower.tail=FALSE)
#     }
#     enr[[jj]]$qenriched.goid[ii] = paste(goids[p.adjust(pp)<.1], collapse = "-")
#     enr[[jj]]$np.enriched[ii] = sum(pp<.01)
#     enr[[jj]]$nq.enriched[ii] = sum(p.adjust(pp)<.1)
#   }
# }
# 
# summary(glm(enr[[1]]$np.enriched ~ clusters$effect + clusters$size))
# summary(glm(enr[[2]]$np.enriched ~ clusters$effect + clusters$size))
# summary(glm(enr[[3]]$np.enriched ~ clusters$effect + clusters$size))
# 
# 
# 
# # co
# print("co")
# tmp = clust.perturb(ints.corum, clustering.algorithm = co.alg, 
#                     noise = 0.1, 
#                     iter = 5)
# tmp$size = clusters$size
# 
# summary(glm(enr[[1]]$np.enriched ~ tmp$reproducibility.J + tmp$size))
# summary(glm(enr[[2]]$np.enriched ~ tmp$reproducibility.J + tmp$size))
# summary(glm(enr[[3]]$np.enriched ~ tmp$reproducibility.J + tmp$size))

