source("functions.R")

# get clusters
fn = "../data/clusters_calcEA.txt"
data = read_tsv(fn)
I = data$iter==1
clusters0 = data[I,]


# get interactomes
fn = "../data/interactomes_calcEA.txt"
interactomes = read_tsv(fn)
unqprots = unique(c(interactomes$protA, interactomes$protB))


# get corum edges
fn = "../data/allComplexes.txt"
corum = read_tsv(fn)
corum.edges = character(10^6)
cc = 0
for (ii in 1:nrow(corum)) {
  prots = unlist(strsplit(corum$`subunits(UniProt IDs)`[ii], ";"))
  for (jj in 1:length(prots)) {
    for (kk in 1:length(prots)) {
      if (jj>=kk) next
      cc = cc+1
      tmp = sort(c(prots[jj], prots[kk]))
      corum.edges[cc] = paste(tmp,collapse="-")
    }
  }
}
corum.edges = corum.edges[1:cc]


#
clusters0$Ar = numeric(nrow(clusters0))
clusters0$size = numeric(nrow(clusters0))
clusters0$corum_score = numeric(nrow(clusters0))
clusters0$Ar.zscore = numeric(nrow(clusters0))
for (ii in 1:nrow(clusters0)) {
  print(ii)
  cluster0 = unlist(strsplit(clusters0$cluster[ii], ";"))
  this.algorithm = clusters0$algorithm[ii]
  this.dataset = clusters0$data_type[ii]
  I = interactomes$iter==1 & interactomes$dataset==this.dataset
  interactome0 = interactomes[I,1:2]
  
  # cluster0 as integers
  cluster0_int = match(cluster0, unqprots)
  
  # cluster0 as paired-strings
  nn = length(cluster0)
  cluster0_pairs = character(nn*(nn-1)/2)
  cc = 0
  for (jj in 1:length(cluster0)) {
    for (kk in 1:length(cluster0)) {
      if (jj>=kk) next
      cc = cc+1
      prots = sort(c(cluster0[jj], cluster0[kk]))
      cluster0_pairs[cc] = paste(prots, collapse="-")
    }
  }
  
  # interactome0 as integers
  interactome0_int = data.frame(protA = match(interactome0$protA, unqprots),
                                protB = match(interactome0$protB, unqprots))
  
  # calculate corum_score
  clusters0$corum_score[ii] = sum(cluster0_pairs %in% corum.edges) / length(cluster0_pairs)
  
  # compare cluster0 to all randomized clusters from the same algorithm
  unqIters = 2:100#sort(unique(data$iter[I]))
  Ar = numeric(length(unqIters))
  Ar_z_iter = 1
  Ar_null = matrix(nrow = length(unqIters), ncol = Ar_z_iter)
  for (jj in 1:length(unqIters)) {
    # cluster assignment reproducibility
    I = data$algorithm %in% this.algorithm & data$noise_mag>0 & data$iter==unqIters[jj]
    Ar[jj] = calcA(cluster0, data$cluster[I])
    
    # calculate null distribution
    for (kk in 1:Ar_z_iter) {
      cluster0.fake = unqprots[sample(length(unqprots), length(cluster0))]
      Ar_null[jj,kk] = calcA(cluster0.fake, data$cluster[I])
    }
  }
  
  clusters0$Ar[ii] = mean(Ar, na.rm=T)
  clusters0$Ar.zscore[ii] = (mean(Ar, na.rm=T) - mean(Ar_null, na.rm=T)) / sd(Ar_null, na.rm=T)
}



x = clusters0$corum_score
y = clusters0$Ar.zscore

test = function(x,y) suppressWarnings(cor.test(x, y, method="spearman")$estimate)
rho = test(x,y)                                             # Test statistic
rho.iter = replicate(10^3, test(x, sample(y, length(y))))   # Simulated permutation distribution

p.out = sum(abs(rho.iter) > rho)    # Count of strict (absolute) exceedances
p.at = sum(abs(rho.iter) == rho)    # Count of equalities, if any
(p.out + p.at /2) / length(p) # Proportion of exceedances: the p-value.
