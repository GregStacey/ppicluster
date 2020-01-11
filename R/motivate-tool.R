
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
  
  # store in broken.mat, type.mat
  for (jj in 1:length(clusts)) {
    cc = cc+1
    prots = unlist(strsplit(clusts[jj], ";"))
    broken$id.clust[cc] = jj
    broken$iter[cc] = iter
    broken$Ji[cc] = calcA(clusts[jj], clusts2)
    broken$n.scrambled[cc] = sum(prots %in% scrambledprots)
    broken$frac.scrambled[cc] = broken$n.scrambled[cc] / length(prots)
  }
  save(clusts,broken, file="../data/broken.Rda")
}
save(clusts,broken, file="../data/broken.Rda")
