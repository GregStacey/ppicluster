# predict real world cluster change

source("functions.R")
source("./clust-perturb-tool/clust-perturb.R")
source("./clust-perturb-tool/functions.R")
source("clusterone_java.R")

# load networks
# # corum
# load("../data/corum_old/all_corum_ints.Rda") # ints.old.corum
# all.corum.ints = list()
# for (ii in 1:5) {
#   all.corum.ints[[ii]] = paste(ints.old.corum[[ii]]$protA, ints.old.corum[[ii]]$protB, sep="-")
# }
# co-fractionation datasets
fn = "../data/interactomes_moredata.txt"
ints.c = as.data.frame(read_tsv(fn))
ints.c = ints.c[ints.c$noise_mag==0, ]
ints.c$ppi = paste(ints.c$protA, ints.c$protB, sep="-")
ints.c = ints.c[ints.c$noise_mag==0,]



# # describe network changes
# # gain and loss between successive CORUM versions
# df = data.frame(edge.gain = numeric(1e3), 
#                         edge.loss = numeric(1e3), 
#                         node.gain = numeric(1e3),
#                         node.loss = numeric(1e3),
#                         n0n0 = numeric(1e3),
#                         n0n1 = numeric(1e3),
#                         n1n1 = numeric(1e3),
#                         nedge1 = numeric(1e3),
#                         nnode1 = numeric(1e3),
#                         network = character(1e3), stringsAsFactors = F)
# for (ii in 2:5) {
#   g1 = graph_from_data_frame(d=ints.old.corum[[ii-1]], directed = F)
#   g2 = graph_from_data_frame(d=ints.old.corum[[ii]], directed = F)
#   df0 = netdiff(g1, g2)
#   df[ii-1,1:7] = df0
#   df$nedge1[ii-1] = length(E(g1))
#   df$nnode1[ii-1] = length(V(g1))
#   df$network[ii-1] = "Successive CORUM releases"
# }
# 
# # gain and loss between replicate pairs
# # load experimental interactomes
# sets = list(c(1:4),c(5:8), 
#             c(9:11),c(12:14), 
#             c(15:17),c(18:20), 
#             c(21:24),c(25:28))
# unqsets = unique(ints.c$dataset)
# cc = 4
# for (ii in 1:length(sets)) {
#   this.set = sets[[ii]]
#   for (jj in 1:length(this.set)) {
#     ia = which(ints.c$dataset==unqsets[this.set[jj]])
#     for (kk in 1:length(this.set)) {
#       if (jj==kk) next
#       cc = cc+1
#       ib = which(ints.c$dataset==unqsets[this.set[kk]])
#       
#       g1 = graph_from_data_frame(d=ints.c[ia,], directed = F)
#       g2 = graph_from_data_frame(d=ints.c[ib,], directed = F)
#       df0 = netdiff(g1, g2)
#       df[cc,1:7] = df0
#       df$nedge1[cc] = length(E(g1))
#       df$nnode1[cc] = length(V(g1))
#       df$network[cc] = "Experiment, rep-rep"
#       if (df$edge.gain[cc]>1e5) error
#     }
#   }
# }
# 
# # gain and loss introduced by experimental noise
# for (ii in 1:length(unqsets)) {
#   ia = ints.c$dataset==unqsets[ii] & ints.c$noise_mag==0
#   ib = ints.c$dataset==unqsets[ii] & ints.c$noise_mag==0.1
#   
#   if (sum(ia)<100 & sum(ib)<100) next
#   cc = cc+1
#   g1 = graph_from_data_frame(d=ints.c[ia,], directed = F)
#   g2 = graph_from_data_frame(d=ints.c[ib,], directed = F)
#   
#   df0 = netdiff(g1, g2)
#   df[cc,1:7] = df0
#   df$nedge1[cc] = length(E(g1))
#   df$nnode1[cc] = length(V(g1))
#   df$network[cc] = "Experiment, noise-induced (10%)"
# }
# df = df[1:cc,]
# 
# # make gain/loss relative to net1
# df$frac.edge.loss = df$edge.loss / df$nedge1
# df$frac.edge.gain = df$edge.gain / df$nedge1
# df$frac.node.loss = df$node.loss / df$nedge1
# df$frac.node.gain = df$node.gain / df$nedge1
# 
# # are new edges preferentially between n2-n2 nodes? n1-n1? n1-n2?
# df$f.old = df$n0n0 / (df$nnode1^2)
# df$f.mix = df$n0n1 / (df$nnode1 * df$node.gain)
# df$f.new = df$n1n1 / (df$node.gain^2)
# 
# ggplot(df, aes(x=edge.loss, y=edge.gain)) + 
#   geom_abline(linetype = "dashed") +
#   facet_wrap(~network, scales = "free") + 
#   geom_point(alpha = .6) + theme_bw() + 
#   xlab("Interactions lost") + ylab("Interactions gained") + ggtitle("edge gain and loss")
# 
# ggplot(df, aes(x=node.loss, y=node.gain)) + 
#   geom_abline(linetype = "dashed") +
#   facet_wrap(~network, scales = "free") + 
#   geom_point(alpha = .6) + theme_bw() + 
#   xlab("Nodes lost") + ylab("Nodes gained") + ggtitle("node gain and loss")
# 
# I = df$frac.edge.gain<3 & df$frac.edge.loss<3
# ggplot(df[I,], aes(x=frac.edge.loss, y=frac.edge.gain)) + 
#   geom_abline(linetype = "dashed") +
#   facet_wrap(~network, scales = "free") + 
#   geom_point(alpha = .6) + theme_bw() + 
#   xlab("Interactions lost") + ylab("Interactions gained") + 
#   ggtitle("fraction edge gain and loss")
# 
# I = df$frac.node.loss<1 & df$frac.node.gain<1
# ggplot(df[I,], aes(x=frac.node.loss, y=frac.node.gain)) + 
#   geom_abline(linetype = "dashed") +
#   facet_wrap(~network, scales = "free") + 
#   geom_point(alpha = .6) + theme_bw() + 
#   xlab("Nodes lost") + ylab("Nodes gained") + 
#   ggtitle("fraction node gain and loss")
# 
# 
# df2 = data.frame(frac = c(df$f.old, df$f.mix, df$f.new),
#                  type = c(rep("n1-n1",nrow(df)), rep("n1-n2",nrow(df)), rep("n2-n2", nrow(df))),
#                  network = rep(df$network, 3), stringsAsFactors = F)
# df2$frac[df2$type=="n1-n1" & grepl("CORUM",df2$network)] = NA
# ggplot(df2, aes(x=frac, fill=type)) + geom_density(alpha=.5) + 
#   facet_wrap(~network, scales="free") +
#   ggtitle("where are new edges in the adjacency matrix? (fraction)")
# 
# df3 = data.frame(count = c(df$n0n0, df$n0n1, df$n1n1),
#                  type = c(rep("n1-n1",nrow(df)), rep("n1-n2",nrow(df)), rep("n2-n2", nrow(df))),
#                  network = rep(df$network, 3), stringsAsFactors = F)
# ggplot(df3, aes(x=log10(count), fill=type)) + geom_density(alpha=.5) + 
#   facet_wrap(~network, scales="free") +
#   ggtitle("where are new edges in the adjacency matrix? (count)")
# 
# ggplot(df, aes(x=(node.gain+node.loss)/nnode1, fill=network)) + 
#   geom_density(alpha = .5) + coord_cartesian(xlim = c(0, 2.5))



## Okay
# after all that, just do replicates!
sets = list(c(1:4),c(5:8), 
            c(9:11),c(12:14), 
            c(15:17),c(18:20), 
            c(21:24),c(25:28))
cc = 0
df.predrep = data.frame(predJ = numeric(1e5), # from clust.perturb
                        repJ = numeric(1e5), # comparing to other set
                        set1 = character(1e5),
                        set2 = character(1e5), stringsAsFactors = F)
for (ii in 1:length(sets)) {
  nn = length(sets[[ii]])
  x1 = sample(sets[[ii]], floor(nn/2))
  x2 = sets[[ii]][!sets[[ii]] %in% x1]
  x1 = sapply(x1, FUN = function(x) paste("us_",x,sep=""))
  x2 = sapply(x2, FUN = function(x) paste("us_",x,sep=""))
  
  # make interactome from subset of replicates
  ia = ints.c$dataset %in% x1 & ints.c$noise_mag==0
  ib = ints.c$dataset %in% x2 & ints.c$noise_mag==0
  protsa = unique(c(ints.c$protA[ia], ints.c$protB[ia]))
  protsb = unique(c(ints.c$protA[ib], ints.c$protB[ib]))
  net1 = ints.c[ia,1:2]
  net2 = ints.c[ib,1:2]
  
  # cluster with clust.perturb
  alg.names = c("k-Med", "MCL", "walktrap", "CO")
  alg = function(x) clusteroneR(x, pp=500, density_threshold = 0.1, java_path = "../java/cluster_one-1.0.jar")
  clust1 = clust.perturb2(net1, clustering.algorithm = alg, noise = 0.15, iters = 5)
  clust2 = clust.perturb2(net2, clustering.algorithm = alg, noise = 0.15, iters = 4)
  
  # calc diff b/w sets
  clust1$repJ.clust2 = numeric(nrow(clust1))
  for (ii in 1:nrow(clust1)) {
    clust1$repJ.clust2[ii] = calcA(clust1$cluster[ii], clust2$cluster)
  }
  clust2$repJ.clust1 = numeric(nrow(clust2))
  for (ii in 1:nrow(clust2)) {
    clust2$repJ.clust1[ii] = calcA(clust2$cluster[ii], clust1$cluster)
  }
  
  # store in dataframe
  I = (cc+1) : (cc+nrow(clust1))
  df.predrep$predJ[I] = clust1$reproducibility.J
  df.predrep$repJ[I] = clust1$repJ.clust2
  df.predrep$set1[I] = paste(x1, collapse = ";")
  df.predrep$set2[I] = paste(x2, collapse = ";")
  cc = cc+nrow(clust1)
  I = (cc+1) : (cc+nrow(clust2))
  df.predrep$predJ[I] = clust2$reproducibility.J
  df.predrep$repJ[I] = clust2$repJ.clust1
  df.predrep$set1[I] = paste(x2, collapse = ";")
  df.predrep$set2[I] = paste(x1, collapse = ";")
  cc = cc+nrow(clust2)
}
df.predrep = df.predrep[1:cc,]
write_tsv(df.predrep, path = "../data/pred_J_expreps.txt")
