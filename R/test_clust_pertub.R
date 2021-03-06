
# load clust.perturb
#source("/Users/gregstacey/Academics/Foster/clust-perturb/R/clust-perturb.R")
#source("/Users/gregstacey/Academics/Foster/clust-perturb/R/functions.R")
source("./clust-perturb-tool/clust-perturb.R")
source("./clust-perturb-tool/functions.R")
source("clusterone_java.R")

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
if (1) {
  load("../data/test_clust_perturb_4algs.Rda")
} else {
  iters = 5
  alg.names = c("k-Med", "MCL", "walktrap", "CO")
  alg = c(function(x) pam(x, 1500),
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
}


# plot results
nn = 1e6
dfrepj = data.frame(repJ = numeric(nn),
                    norm.rank = numeric(nn),
                    noise_mag = numeric(nn),
                    algorithm = character(nn), stringsAsFactors = F)
cc = 0
for (ii in 1:length(noise.range)) {
  I = (cc+1):(cc+nrow(clusters.kmed[[ii]]))
  dfrepj$repJ[I] = clusters.kmed[[ii]]$reproducibility.J
  dfrepj$norm.rank[I] = rank(clusters.kmed[[ii]]$reproducibility.J) / length(I)
  dfrepj$noise_mag[I] = noise.range[ii]
  dfrepj$algorithm[I] = "k-Med"
  cc = cc+length(I)
  
  I = (cc+1):(cc+nrow(clusters.mcl[[ii]]))
  dfrepj$repJ[I] = clusters.mcl[[ii]]$reproducibility.J
  dfrepj$norm.rank[I] = rank(clusters.mcl[[ii]]$reproducibility.J) / length(I)
  dfrepj$noise_mag[I] = noise.range[ii]
  dfrepj$algorithm[I] = "MCL"
  cc = cc+length(I)
  
  I = (cc+1):(cc+nrow(clusters.walk[[ii]]))
  dfrepj$repJ[I] = clusters.walk[[ii]]$reproducibility.J
  dfrepj$norm.rank[I] = rank(clusters.walk[[ii]]$reproducibility.J) / length(I)
  dfrepj$noise_mag[I] = noise.range[ii]
  dfrepj$algorithm[I] = "walktrap"
  cc = cc+length(I)
  
  I = (cc+1):(cc+nrow(clusters.co[[ii]]))
  dfrepj$repJ[I] = clusters.co[[ii]]$reproducibility.J
  dfrepj$norm.rank[I] = rank(clusters.co[[ii]]$reproducibility.J) / length(I)
  dfrepj$noise_mag[I] = noise.range[ii]
  dfrepj$algorithm[I] = "CO"
  cc = cc+length(I)
}
dfrepj = dfrepj[1:cc, ]

ggplot(dfrepj, aes(x=norm.rank, y=repJ, color=as.factor(noise_mag))) +
  facet_wrap(~algorithm) + geom_line()



### how many iterations do we need?



### how does repJ vary across iterations? as a function of noise?
# (hint: it's good if it's consistent!!)



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





### motivate the tool
# show there are clusters that tend to break

load("../data/broken.Rda") # broken

# add cluster size to broken
clust.size = sapply((sapply(clusts, strsplit, ";")), length)
broken$size = clust.size[broken$id.clust]

# look for clusters that, when segregated, tend to be broken
I = broken$n.scrambled==0
ggplot(broken[I,], aes(x=id.clust,y=Ji)) + geom_point(alpha=.2)
I = broken$frac.scrambled<=0.05
ggplot(broken[I,], aes(x=id.clust,y=Ji)) + geom_jitter(alpha=.2, height = 0.02)

# look for examples
id = unique(broken$id.clust[broken$size==10])
I = broken$n.scrambled==0 & broken$id.clust%in%id
ggplot(broken[I,], aes(x=as.character(id.clust),y=Ji)) + geom_violin(alpha=.2)

# make this plot:
#   x = cluster rank
#   y = Ji
#   color = n.broken
#   facet = size
df = data.frame(Ji = numeric(1e5), rank = numeric(1e5), norm.rank = numeric(1e5),
                size = numeric(1e5), n.scrambled = numeric(1e5),
                stringsAsFactors = F)
cc = 0
#unqsize = c(5, 10, 25, 50)
unqsize = c(1e5)
scramblerange = c(0, 1, 2)
for (ii in 1:length(unqsize)) {
  for (jj in 1:length(scramblerange)) {
    I = broken$size<=unqsize[ii] & broken$n.scrambled==scramblerange[jj]
    if (jj==3) I = broken$size<=unqsize[ii] & broken$n.scrambled>=scramblerange[jj]
    if (sum(I)==0) next
    I2 = (cc+1) : (cc+sum(I))
    df$Ji[I2] = sort(broken$Ji[I], decreasing = F)
    df$rank[I2] = 1:sum(I)
    df$norm.rank[I2] = seq(from=0, to=1, length=sum(I))
    df$size[I2] = unqsize[ii]
    df$n.scrambled[I2] = scramblerange[jj]
    cc = cc+sum(I)
  }
}
df = df[1:cc,]
ggplot(df, aes(x=norm.rank, y=Ji, color=as.character(n.scrambled))) + 
  facet_wrap(~size) + geom_line(alpha=.99)

ggplot(df, aes(y=Ji, x=as.character(n.scrambled))) + 
  facet_wrap(~size) + geom_violin(alpha=.99)

# Conclusion: Yes! There are clusters that, when segregated from noise, still tend to break.




### AHA!
# the tool should produce a connection matrix, where each edge is measuring
# how likely that pair was to show up in noised clusters.
# do this for every original cluster.
# use it to score clusters and show if/how they break apart.
# identify a cluster that splits in half sensibly.
# conversely, show a well-connected cluster that survives many edge shufflings.
#
# questions
# 1. how to include this noise-segregating/noise-adding thing to the tool?
# 2. how to add noise in the tool? why not randomly choose n.add and n.remove?
#    i.e. not require strict shuffling
#    or do you want to directly probe each cluster?
#
# by segregating noise this way, you get 2 values per clust: n.scrambled + Ji.
# you should be able to fit a model this data to see if it's higher/lower than random!





### predict release-to-release variability of corum clusters

# load corums
load("../data/corum_old/all_corum_ints.Rda")

# calculate reproducibility of each release
if (1) {
  load("../data/corum_old/clusters_tool.Rda") # clusters.old.corum
} else {
  clusters.old.corum = list()
  co.alg = function(x) neR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar")
  for (ii in 1:length(ints.old.corum)) {
    clusters.old.corum[[ii]]  = clust.perturb(ints.old.corum[[ii]], 
                                              clustering.algorithm = co.alg, 
                                              noise = 0.1, 
                                              iter = 1,
                                              edge.list.format = NULL, 
                                              cluster.format = NULL)
  }
  save(clusters.old.corum, file = "../data/corum_old/clusters_tool.Rda")
}

# calculate release-to-release variability
for (ii in 1:(length(clusters.old.corum)-1)) {
  print(ii)
  clusters.old.corum[[ii]]$prev.release.J = rep(NA, nrow(clusters.old.corum[[ii]]))
  these.clusters = clusters.old.corum[[ii]]$cluster
  ref.clusters = clusters.old.corum[[ii+4]]$cluster
  for (jj in 1:length(these.clusters)) {
    clusters.old.corum[[ii]]$prev.release.J[jj] = calcA(these.clusters[jj], ref.clusters)
  }
  
  plot(clusters.old.corum[[ii]]$prev.release.J, 
       clusters.old.corum[[ii]]$reproducibility.J)
}



# ... okay, this doesn't work great.
#
# works well for clusters that aren't reproduced in the next corum,
# but doesn't predict clusters that ARE reproduced.
#
# that is, the tool is too pessimistic about reproducibility.

# explore...
tmp = clusters.old.corum[[5]]
df = data.frame(bins = seq(from=0, to=.9, length=10),
                nextJ = numeric(10))
for (ii in 1:nrow(df)) {
  ia = tmp$reproducibility.J>df$bins[ii] & tmp$reproducibility.J<(df$bins[ii]+0.1)
  df$nextJ[ii] = mean(tmp$prev.release.J[ia])
}
ggplot(df, aes(x=bins, y=nextJ)) + geom_line()

tmp = clusters.old.corum[[1]]
tmp$rep.fact = tmp$reproducibility.J > median(tmp$reproducibility.J)
ggplot(tmp, aes(x=prev.release.J, fill = rep.fact)) + geom_density(alpha=.6)



# okay rewind.
# this only works if one random iteration is predictive of another.
# test that
tmp1 = clust.perturb(ints.old.corum[[ii]], 
                     clustering.algorithm = co.alg, 
                     noise = 0.1, 
                     iter = 1)
tmp2 = clust.perturb(ints.old.corum[[ii]], 
                     clustering.algorithm = co.alg, 
                     noise = 0.1, 
                     iter = 1)
tmp3 = clust.perturb(ints.old.corum[[ii]], 
                     clustering.algorithm = co.alg, 
                     noise = 0.1, 
                     iter = 1)
tmp4 = clust.perturb(ints.old.corum[[ii]], 
                     clustering.algorithm = co.alg, 
                     noise = 0.01, 
                     iter = 1)
cor.test(tmp1$reproducibility.J, tmp2$reproducibility.J) # correlated
cor.test(tmp1$reproducibility.J, tmp3$reproducibility.J) # correlated
cor.test(tmp2$reproducibility.J, tmp3$reproducibility.J) # correlated
cor.test(tmp1$reproducibility.J, tmp4$reproducibility.J) #


# scratchpad
# try to get tool to be predictive of corum breakage

# now.
# if one random iteration is predictive of another, treat corum->corum as an iteration.

# define: how much noise from corum1->5
ints1 = paste(ints.old.corum[[1]]$protA, ints.old.corum[[1]]$protB, sep="-")
ints5 = paste(ints.old.corum[[5]]$protA, ints.old.corum[[5]]$protB, sep="-")
n.added = sum(!ints5 %in% ints1)
n.removed = sum(!ints1 %in% ints5)
print(paste("fraction added =",n.added/length(ints1)))
print(paste("fraction removed =",n.removed/length(ints1)))

# calculate J(corum1->5)
these.clusters = clusters.old.corum[[1]]$cluster
ref.clusters = clusters.old.corum[[5]]$cluster
J15 = numeric(nrow(clusters.old.corum[[1]]))
for (jj in 1:length(these.clusters)) {
  J15[jj] = calcA(these.clusters[jj], ref.clusters)
}

# now, simulate this with random noise
tmp.add = addcorum(ints.old.corum[[1]], n.added/length(ints1))
tmp.remove = removecorum(ints.old.corum[[1]], n.removed/length(ints1))
# find added + removed interactions
i.add = !paste(tmp.add$protA, tmp.add$protB, sep="-") %in% 
  paste(ints.old.corum[[1]]$protA, ints.old.corum[[1]]$protB, sep="-")
i.remove = !paste(ints.old.corum[[1]]$protA, ints.old.corum[[1]]$protB, sep="-") %in% 
  paste(tmp.remove$protA, tmp.remove$protB, sep="-")
# put it together
ints1.noised = rbind(ints.old.corum[[1]][!i.remove,], tmp.add[i.add,])

# check: how much noise from corum1->noise
ints1 = paste(ints.old.corum[[1]]$protA, ints.old.corum[[1]]$protB, sep="-")
intsnoise = paste(ints1.noised$protA, ints1.noised$protB, sep="-")
n.added = sum(!intsnoise %in% ints1)
n.removed = sum(!ints1 %in% intsnoise)
print(paste("fraction added =",n.added/length(ints1)))
print(paste("fraction removed =",n.removed/length(ints1)))

# cluster
clust1 = co.alg(ints.old.corum[[1]])
clust1.noised = co.alg(ints1.noised)
print(paste("are clust1 and clusters.old.corum[[1]] identical?",
            identical(clust1, as.list(clusters.old.corum[[1]]$cluster))))

# calculate J(corum1->noise)
these.clusters = unlist(clust1)
ref.clusters = unlist(clust1.noised)
J1noise = numeric(nrow(clusters.old.corum[[1]]))
for (jj in 1:length(these.clusters)) {
  J1noise[jj] = calcA(these.clusters[jj], ref.clusters)
}

plot(J15, J1noise)

# aha!
# it's the same pattern
# J15 has way more 1 values, and these are throwing it off.
# J15 is positive-skewed, J1noise is negative-skewed.
# that is, there is something inherently different between J1->5 and J1->noise processes.
#
# therefore
# look for differences between (ints1->5)vs(ints1->noise) and (ints5)vs(intsnoise)

# diff1: node loss
# ints1->5 adds 1926 and loses 282 nodes, ints1->noise adds and loses none
nodes1 = unique(c(ints.old.corum[[1]]$protA,ints.old.corum[[1]]$protB))
nodes5 = unique(c(ints.old.corum[[5]]$protA,ints.old.corum[[5]]$protB))
nodesnoise = unique(c(ints1.noised$protA,ints1.noised$protB))
n.gain = sum(!nodes5 %in% nodes1)
n.loss = sum(!nodes1 %in% nodes5)

# oooooooh
# maybe the edges added in ints1->5 are between NEW nodes!
# therefore they don't screw up existing clusters.
# yeah, that makes sense!

# therefore, clust.perturb will be poorly predictive when the primary change is node loss+gain.
# therefore, think of a situation where edges change but nodes stay the same.

# TEST - filter ints and ints5 to common nodes
#      - compare J(corum1.filtered->corum5.filtered) to J(corum1.filtered->corum1.filtered.noised)
# calculate J(corum1.filtered->corum5.filtered)
nodes1 = unique(c(ints.old.corum[[1]]$protA,ints.old.corum[[1]]$protB))
nodes5 = unique(c(ints.old.corum[[5]]$protA,ints.old.corum[[5]]$protB))
nodes = intersect(nodes1, nodes5)
i1 = ints.old.corum[[1]]$protA %in% nodes & ints.old.corum[[1]]$protB %in% nodes
ints1.old.filtered = ints.old.corum[[1]][i1,]
ints1.filtered = paste(ints1.old.filtered$protA, ints1.old.filtered$protB, sep="-")
i5 = ints.old.corum[[5]]$protA %in% nodes & ints.old.corum[[5]]$protB %in% nodes
ints5.old.filtered = ints.old.corum[[5]][i5,]
ints5.filtered = paste(ints5.old.filtered$protA, ints5.old.filtered$protB, sep="-")
clust1.filtered = co.alg(ints1.old.filtered)
clust5.filtered = co.alg(ints5.old.filtered)
# calculate J(corum1.filtered->corum5.filtered)
these.clusters = unlist(clust1.filtered)
ref.clusters = unlist(clust5.filtered)
J1filtered = numeric(length(these.clusters))
for (jj in 1:length(these.clusters)) {
  J1filtered[jj] = calcA(these.clusters[jj], ref.clusters)
}
# define: how much noise from corum1filtered->5filtered
n.added = sum(!ints5.filtered %in% ints1.filtered)
n.removed = sum(!ints1.filtered %in% ints5.filtered)
print(paste("fraction added =",n.added/length(ints1)))
print(paste("fraction removed =",n.removed/length(ints1)))
# simulate the same level of noise (random) 
tmp.add = addcorum(ints1.old.filtered, n.added/sum(i1))
tmp.remove = removecorum(ints1.old.filtered, n.removed/sum(i1))
# find added + removed interactions
i.add = !paste(tmp.add$protA, tmp.add$protB, sep="-") %in% 
  paste(ints1.old.filtered$protA, ints1.old.filtered$protB, sep="-")
i.remove = !paste(ints1.old.filtered$protA, ints1.old.filtered$protB, sep="-") %in% 
  paste(tmp.remove$protA, tmp.remove$protB, sep="-")
# put it together
ints1.filtered.noised = rbind(ints1.old.filtered[!i.remove,], tmp.add[i.add,])
# cluster
clust1.filtered = co.alg(ints1.old.filtered)
clust1.filtered.noised = co.alg(ints1.filtered.noised)
# calculate J(corum1.filtered->corum1.filtered.noised)
these.clusters = unlist(clust1.filtered)
ref.clusters = unlist(clust1.noised)
J1filtered.noise = numeric(length(these.clusters))
for (jj in 1:length(these.clusters)) {
  J1filtered.noise[jj] = calcA(these.clusters[jj], ref.clusters)
}
cor.test(J1filtered, J1filtered.noise)


# YES! It works!!
#
# now do it with the tool
co.alg = function(x) clusteroneR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar")
clusters1.old.corum.filtered  = clust.perturb(ints1.old.filtered, 
                                              clustering.algorithm = co.alg, 
                                              noise = 0.01, 
                                              iter = 15,
                                              edge.list.format = NULL, 
                                              cluster.format = NULL)
print(paste("are clust1.filtered and clusters1.old.corum.filtered identical?",
            identical(clust1.filtered, as.list(clusters1.old.corum.filtered$cluster))))
cor.test(J1filtered, clusters1.old.corum.filtered$reproducibility.J)


# Nooooo...
# it failed...
#
# try with separate add+remove levels
# now do it with the tool
co.alg = function(x) clusteroneR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar")
clusters1.old.corum.filtered2  = clust.perturb2(ints1.old.filtered, 
                                                clustering.algorithm = co.alg, 
                                                add.noise = n.added/sum(i1), 
                                                remove.noise = n.removed/sum(i1), 
                                                iter = 3,
                                                edge.list.format = NULL, 
                                                cluster.format = NULL)
print(paste("are clust1.filtered and clusters1.old.corum.filtered identical?",
            identical(clust1.filtered, as.list(clusters1.old.corum.filtered2$cluster))))
cor.test(J1filtered, clusters1.old.corum.filtered2$reproducibility.J)






##### Waaaaait.
# Why is low noise repJ so different from segJ for size=3 clusters?
# In theory they should be kind of the same at low noise!
#   from figure05.R
#     mean(jmat[df$size==3,]) = 0.42
#     mean(df$seg0.repJ[df$size==3], na.rm=T) = 0.99
co.alg = function(x) clusteroneR(x, pp=500, density_threshold = 0.1, 
                                 java_path = "../java/cluster_one-1.0.jar")
test.clust = clust.perturb2(ints.corum, clustering.algorithm = co.alg, noise=0.1, iters=1)





##### okay wtheck.
#   - clust.perturb is not reproducing mcl and kmed results from fig3
#   - in particular, Ji is way too low
#   - also, clusters.kmed$cluster has NAs in it
#   - also, clusters from non-noised network don't match

unqprots = unique(c(ints.corum$protA, ints.corum$protB))

alg = function(x) pam(x, 100)
edge.list.format = pam.edge.list.format

# start with simple network + first part of tool
ints = ints.corum[1:1000, ]
unqprots = unique(c(ints$protA, ints$protB))
cluster.format = function(x) pam.cluster.format(x,unqprots = unqprots)
network.input = edge.list.format(ints)
tmp = clustering.algorithm(network.input)
ctool.part = unlist(cluster.format(tmp))
cfig3.part = pamclust(ints, 100)
identical(ctool.part, cfig3.part)
# AHA!!!!! 
# it's the unqprots. it must be for the current network!

# next, simple network + entire tool
ints = ints.corum[1:1000, ]
cluster.format = function(x) pam.cluster.format(x,unqprots = unqprots)
unqprots = unique(c(ints$protA, ints$protB))
ctool.part = clust.perturb2(ints, clustering.algorithm = alg, noise = 0.1, iter = 2, 
                            edge.list.format = pam.edge.list.format, cluster.format = cluster.format)
cfig3.part = pamclust(ints, 100)
identical(ctool.part$cluster, cfig3.part)
# WOOHOO!
# Yup, it's the unqprots.

# now do all of corum
ints = ints.corum
alg = function(x) pam(x, 1500)
unqprots = unique(c(ints$protA, ints$protB))
cluster.format = function(x) pam.cluster.format(x,unqprots = unqprots)
ctool.all = clust.perturb2(ints, clustering.algorithm = alg, noise = 0.1, iter = 2, 
                           edge.list.format = pam.edge.list.format, cluster.format = cluster.format)
cfig3.all = pamclust(ints, 1500)
identical(ctool.all$cluster, cfig3.all)
# YAY!
# Still identical! Okay, I fixed that!!

# Now, why is Ji so low?
# try both ways with a lower noise (5%)
ints = ints.corum
noise = 0.05
# tool
alg = function(x) pam(x, 1500)
unqprots = unique(c(ints$protA, ints$protB))
cluster.format = function(x) pam.cluster.format(x,unqprots = unqprots)
ctool.all = clust.perturb2(ints, clustering.algorithm = alg, noise = noise, iter = 1, 
                           edge.list.format = pam.edge.list.format, cluster.format = cluster.format)
# fig 3 method
tmp = data.frame(cluster = character(1e5),
                 noise = numeric(1e5),
                 Ji = numeric(1e5), stringsAsFactors = F)
clustpam0 = pamclust(ints.corum, nclust = 1500)
ints.shuffle = shufflecorum(ints.corum, noise)
clustpam1 = pamclust(ints.shuffle, nclust = 1500)
Jfig3 = numeric(length(clustpam0))
for (ii in 1:length(Jfig3)) {
  tmp = calcA2(clustpam0[ii], clustpam1)
  Jfig3[ii] = tmp[[1]]
}
identical(clustpam0, ctool.all$cluster)
t.test(Jfig3, ctool.all$reproducibility.J)
cor.test(Jfig3, ctool.all$reproducibility.J)
# OOF
# yeah, the clusters are identical, but
# Ji is very different between the two methods

# is edge.list.format() screwing up the network?
edge.list.format = pam.edge.list.format
e0 = edge.list.format(ints.corum[1:1000, ])
e1 = edge.list.format(shufflecorum(ints.corum[1:1000, ], 0.05))

# return noise.clusts from the tool and manually calculate Ji
ints = ints.corum
noise = 0.05
# tool
alg = function(x) pam(x, 1500)
unqprots = unique(c(ints$protA, ints$protB))
cluster.format = function(x,y) pam.cluster.format(x,y = unqprots)
tmp = clust.perturb2(ints, clustering.algorithm = alg, noise = noise, iter = 1, 
                     edge.list.format = pam.edge.list.format, cluster.format = cluster.format)
# ints.shuffle is different...
e0 = edge.list.format(ints.corum)
sum(tmp[[3]]==0 & e0==0) / nrow(ints.corum)


# OHHHHHH
# maybe shuffling corum can remove a node or two
# even losing 1 would change unqprots
# try again when recalculating unqprots after shuffling
ints = ints.corum[1:2000,]
noise = 0.05
# tool
alg = function(x) pam(x, 100)
unqprots = unique(c(ints$protA, ints$protB))
cluster.format = function(x) pam.cluster.format(x,unqprots = unqprots)
ctool.all = clust.perturb2(ints, clustering.algorithm = alg, noise = noise, iter = 1, 
                           edge.list.format = pam.edge.list.format, cluster.format = cluster.format)
# fig 3 method
tmp = data.frame(cluster = character(1e5),
                 noise = numeric(1e5),
                 Ji = numeric(1e5), stringsAsFactors = F)
clustpam0 = pamclust(ints, nclust = 100)
ints.shuffle = shufflecorum(ints, noise)
clustpam1 = pamclust(ints.shuffle, nclust = 100)
Jfig3 = numeric(length(clustpam0))
for (ii in 1:length(Jfig3)) {
  tmp = calcA2(clustpam0[ii], clustpam1)
  Jfig3[ii] = tmp[[1]]
}
# recalculate J from tool
ints.shuffle = shufflecorum(network, noise)
unqprots = unique(c(ints.shuffle$protA, ints.shuffle$protB))
ints.shuffle = edge.list.format(ints.shuffle)
these.clusters = clustering.algorithm(ints.shuffle)
these.clusters = unlist(pam.cluster.format(these.clusters, unqprots = unqprots))
Jrecalc = numeric(length(clustpam0))
for (ii in 1:length(clustpam0)) {
  Jrecalc[ii] = calcA(clustpam0[ii], these.clusters)
}
identical(clustpam0, ctool.all$cluster)
t.test(Jfig3, Jrecalc)
cor.test(Jfig3, Jrecalc)


# okay
# I replaced cluster.format() with pam.cluster.format()
# and recalculated unqprots after shuffling. Let's do this thang!
ints = ints.corum[1:10000,]
noise = 0.1
# tool
alg = function(x) pam(x, 250)
unqprots = unique(c(ints$protA, ints$protB))
cluster.format = function(x, y) pam.cluster.format(clusts = x,unqprots = y)
ctool.all = clust.perturb2(ints, clustering.algorithm = alg, noise = noise, iter = 5, 
                           edge.list.format = pam.edge.list.format, cluster.format = cluster.format)
# fig 3 method
tmp = data.frame(cluster = character(1e5),
                 noise = numeric(1e5),
                 Ji = numeric(1e5), stringsAsFactors = F)
clustpam0 = pamclust(ints, nclust = 250)
ints.shuffle = shufflecorum(ints, noise)
clustpam1 = pamclust(ints.shuffle, nclust = 250)
Jfig3 = numeric(length(clustpam0))
for (ii in 1:length(Jfig3)) {
  tmp = calcA2(clustpam0[ii], clustpam1)
  Jfig3[ii] = tmp[[1]]
}
identical(clustpam0, ctool.all$cluster)
t.test(Jfig3, ctool.all$reproducibility.J)
cor.test(Jfig3, ctool.all$reproducibility.J)






# oi vey.
# do the same process ^ for mcl
# except you can't fully replicate Figure 3, since that data's from matlab mcl
ints = ints.corum[1:2000,]
noise = 0.01
# tool
alg = function(x) mcl(x, addLoops = FALSE)
unqprots = unique(c(ints$protA, ints$protB))
#cluster.format = function(x) mcl.cluster.format(x,unqprots = unqprots)
ctool.all = clust.perturb2(ints, clustering.algorithm = alg, noise = noise, iter = 1, 
                           edge.list.format = mcl.edge.list.format, cluster.format = mcl.cluster.format)
# fig 3 method
tmp = data.frame(cluster = character(1e5),
                 noise = numeric(1e5),
                 Ji = numeric(1e5), stringsAsFactors = F)
clustmcl0 = unlist(mcl.cluster.format(mcl(mcl.edge.list.format(ints), addLoops = F)))
ints.shuffle = shufflecorum(ints, noise)
clustmcl1 = unlist(mcl.cluster.format(mcl(mcl.edge.list.format(ints.shuffle), addLoops = F)))
Jfig3 = numeric(length(clustmcl0))
for (ii in 1:length(Jfig3)) {
  tmp = calcA2(clustmcl0[ii], clustmcl1)
  Jfig3[ii] = tmp[[1]]
}
identical(clustmcl0, ctool.all$cluster)
t.test(Jfig3, ctool.all$reproducibility.J)
cor.test(Jfig3, ctool.all$reproducibility.J)
# con

# 1. is shuffling working
ints = ints.corum[1:2000,]
noise = 0.1
ints.shuffle = shufflecorum(ints, noise)
print(sum(ints$protA==ints.shuffle$protA & ints$protB==ints.shuffle$protB) / nrow(ints))

# 2. is edge formatting good?
edges0 = mcl.edge.list.format(ints)
edges1 = mcl.edge.list.format(ints.shuffle)
sum(rownames(edges0) == rownames(edges1)) / length(rownames(edges0))
# AHA!
# the rownames are getting shuffled
# which means the protein order is changing
# try with passing unqprots...
unqprots0 = unique(c(ints$protA, ints$protB))
unqprots1 = unique(c(ints.shuffle$protA, ints.shuffle$protB))
edges0 = mcl.edge.list.format(ints, unqprots = unqprots0)
edges1 = mcl.edge.list.format(ints.shuffle, unqprots = unqprots1)
sum(rownames(edges0) == rownames(edges1)) / length(rownames(edges0))
# break it up
G0 = graph.data.frame(ints,directed=FALSE)
A0 = as_adjacency_matrix(G0,type="both",names=TRUE,sparse=FALSE)
G1 = graph.data.frame(ints.shuffle,directed=FALSE)
A1 = as_adjacency_matrix(G1,type="both",names=TRUE,sparse=FALSE)
sum(rownames(A0) == rownames(A1)) / length(rownames(A0))
# yeah, proteins are getting reordered
# and I don't understand what G is
# try forcing order of A
unqprots0 = unique(c(ints$protA, ints$protB))
G0 = graph.data.frame(ints,directed=FALSE)
A0 = as_adjacency_matrix(G0,type="both",names=TRUE,sparse=FALSE)
I0 = match(rownames(A0), unqprots0)
G1 = graph.data.frame(ints.shuffle,directed=FALSE)
A1 = as_adjacency_matrix(G1,type="both",names=TRUE,sparse=FALSE)
I1 = match(unqprots0, rownames(A1))
sum(rownames(A0) == rownames(A1[I1, I1])) / length(rownames(A0))
# Yes! that works
# now try putting that in the tool
ints = ints.corum[1:4000,]
noise = 0.2
ctool.all = clust.perturb2(ints, clustering.algorithm = alg, 
                           noise = noise, iter = 1, 
                           edge.list.format = mcl.edge.list.format, 
                           cluster.format = mcl.cluster.format)

# wait... it doesn't work
# why
ints = ints.corum0[10001:21000,]
ints.shuffle = shufflecorum(ints, 0.1)
unqprots0 = unique(c(ints$protA, ints$protB))
G = graph.data.frame(ints.shuffle,directed=FALSE)
A = as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)
# A doesn't have all the edges. why?
for (ii in 1:nrow(ints.shuffle)) {
  ia = which(rownames(A)==ints.shuffle$protA[ii])
  ib = which(colnames(A)==ints.shuffle$protB[ii])
  if (A[ia,ib]<1) error
}
# correct protein ordering
I = match(unqprots0, rownames(A))
A = A[I,I]
A[is.na(A)] = 0
rownames(A) = colnames(A) = unqprots0


#...
# okay, try it now with the new code

# 1. is shuffling good?
ints = ints.corum0[1:2000,]
noise = 0.1
ints.shuffle = shufflecorum(ints, noise)
print(sum(ints$protA==ints.shuffle$protA & ints$protB==ints.shuffle$protB) / nrow(ints))

# 2. is edge formatting good?
x0 = mcl.edge.list.format(ints)
x1 = mcl.edge.list.format(ints.shuffle)

# 3. is clustering good
y0 = mcl(x0,addLoops = F)
y1 = mcl(x1,addLoops = F)

# 4. cluster formatting is good
z0 = mcl.cluster.format(y0)
z1 = mcl.cluster.format(y1)

# 5. does everything work?
alg = function(x) mcl(x, addLoops = F)
noise = 0.25
ctool = clust.perturb2(ints, clustering.algorithm = alg, 
                       noise = noise, iter = 3, 
                       edge.list.format = mcl.edge.list.format, 
                       cluster.format = mcl.cluster.format)


# argh
# things look okay
# but mcl is failing lots!
# see if another package would do better...
require(hbm)
ints = ints.corum[1:5000,]
ints1 = shufflecorum(ints, 0.1)
adj = mcl.edge.list.format(ints)
adj1 = mcl.edge.list.format(ints1)
clust = hbm::mcl(adj, 2)
clust1 = hbm::mcl(adj1, 2)
# test it with noise
ints = ints.corum[1:5000,]
ints1 = shufflecorum(ints, 0.2)
unqprots = unique(c(ints$protA, ints$protB))
unqprots1 = unique(c(ints1$protA, ints1$protB))
adj = mcl.edge.list.format(ints,unqprots)
adj1 = mcl.edge.list.format(ints1)
clust = mcl.cluster.format(hbm::mcl(adj, 2), unqprots)
clust1 = mcl.cluster.format(hbm::mcl(adj1, 2), unqprots1)
JJ = numeric(length(clust))
for (ii in 1:length(clust)) {
  JJ[ii] = calcA(clust[ii], clust1)
}
# check the tool works
alg = function(x) mymcl(x, infl = 2)
ctool = clust.perturb2(ints, clustering.algorithm = alg, 
                       noise = noise, iters = 3, 
                       edge.list.format = mcl.edge.list.format, 
                       cluster.format = mcl.cluster.format)

# argh, try it again with my change to hbm::mcl
alg = function(x) mymcl(x, infl = 2, iter = 20, verbose = T)
noise = 0.1
iters = 3
#
ints = ints.corum[1:1000,]
unqprots = unique(c(ints$protA, ints$protB))
x0 = mcl.edge.list.format(ints)
y0 = alg(x0)
z0 = mcl.cluster.format(y0, unqprots = unqprots)
#
ints1 = shufflecorum(ints, 0.1)
unqprots1 = unique(c(ints1$protA, ints1$protB))
x1 = mcl.edge.list.format(ints1)
y1 = alg(x1)
z1 = mcl.cluster.format(y1, unqprots = unqprots1)
JJ = numeric(length(z0))
for (ii in 1:length(JJ)) {
  JJ[ii] = calcA(z0[ii], z1)
}
# try with the tool
ctool = clust.perturb2(ints, clustering.algorithm = alg, 
                       noise = noise, iters = 3, 
                       edge.list.format = mcl.edge.list.format, 
                       cluster.format = mcl.cluster.format)
# try a higher noise range...
ctool = clust.perturb2(ints, clustering.algorithm = alg, 
                       noise = 0.3, iters = 3, 
                       edge.list.format = mcl.edge.list.format, 
                       cluster.format = mcl.cluster.format)
# try all of corum with the tool...
alg = function(x) mymcl(x, infl = 2, iter = 2, verbose = T)
ctool = clust.perturb2(ints.corum, clustering.algorithm = alg, 
                       noise = 0.1, iters = 1, 
                       edge.list.format = mcl.edge.list.format, 
                       cluster.format = mcl.cluster.format)


