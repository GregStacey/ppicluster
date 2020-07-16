
# install packages if necessary
list.of.packages = c("tidyverse", "ggplot2", "magrittr", "readr", "cluster", "tidyr", "NMI",
                     "fossil", "infotheo", "igraph", "dils", "hbm", "tools", "flavin", 
                     "resolution", "leiden", "ProNet", "dbscan", "clustperturb")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

require(tidyverse)
require(ggplot2)
require(magrittr)
require(readr)
require(reshape2)
require(viridis)
require(cluster)
require(tidyr)
require(dplyr)
require(NMI)
require(fossil)
#require(VennDiagram)
require(infotheo)
require(igraph)
require(dils)
#require(MCL)
require(hbm)
require(tools)
mcl = MCL::mcl
#require(ontologyIndex)
require(flavin)
require(resolution) # louvain clustering with resolution
require(leiden)
require(ProNet)
require(dbscan)
require(clust.perturb)
source("clusterone_java.R")
source("mymcl.R")


blank_theme = theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(), legend.position="none")

# define functions
calcA = function(this.cluster, clusters) {
  # J, Ji
  # assignment reproducibility
  # essentially the maximum Jaccard index
  #
  # this.cluster = single string, with semicolon-delimited IDs
  # clusters = vector of strings, all semicolon-delimited IDs
  
  this.cluster = unlist(strsplit(this.cluster, ";"))
  JJ = numeric(length(clusters))
  for (ii in 1:length(clusters)) {
    that.cluster = unlist(strsplit(clusters[ii], ";"))
    JJ[ii] = length(intersect(this.cluster, that.cluster)) / 
      length(unique(c(this.cluster, that.cluster)))
  }
  A = max(JJ, na.rm=T)
  return(A)
}


calcAz = function(this.cluster, clusters) {
  # z-score of J, Ji
  this.cluster = unlist(strsplit(this.cluster, ";"))
  JJ = numeric(length(clusters))
  for (ii in 1:length(clusters)) {
    that.cluster = unlist(strsplit(clusters[ii], ";"))
    JJ[ii] = length(intersect(this.cluster, that.cluster)) / 
      length(unique(c(this.cluster, that.cluster)))
  }
  A = max(JJ, na.rm=T)
  
  # now do this 1000 more times with random clusters
  # convert clusters to a data.frame = [node, cluster.id]
  iterMax = 50
  df.cluster = clusters.to.df(clusters)
  A.rand = rep(NA, iterMax)
  for (iter in 1:iterMax) {
    #print(iter)
    JJ.rand = numeric(length(clusters))
    rand.cluster = df.cluster
    rand.cluster$cluster = rand.cluster$cluster[sample(nrow(df.cluster), nrow(df.cluster))]
    for (ii in 1:length(clusters)) {
      that.cluster = rand.cluster$prots[rand.cluster$cluster==ii]
      JJ.rand[ii] = length(intersect(this.cluster, that.cluster)) / 
        length(unique(c(this.cluster, that.cluster)))
    }
    A.rand[iter] = max(JJ.rand, na.rm=T)
  }
  
  Az = (A - mean(A.rand)) / sd(A.rand)
  tmp = list(Az, mean(A.rand), sd(A.rand))
  return(tmp)
}


clusters.to.df = function(clusters) {
  prots = unlist(strsplit(clusters, ";"))
  df = data.frame(prots = prots, stringsAsFactors = F)
  df$cluster = numeric(nrow(df))
  cc = 0
  for (ii in 1:length(clusters)) {
    nn = length(unlist(strsplit(clusters[ii], ";")))
    I = (cc+1) : (cc+nn)
    df$cluster[I] = ii
    cc = cc+nn
  }
  return(df)
}

df.to.cluster = function(df) {
  unqid = unique(df$cluster)
  clusters = character(length(unqid))
  for (ii in 1:length(unqid)) {
    clusters[ii] = paste(unique(df$prots[df$cluster == unqid[ii]]), collapse = ";")
  }
  return(clusters)
}


calcE = function(cluster0_int, interactome0_int, this.interactome) {
  # edge reproducibility
  #
  # cluster0 = integer vector (cluster of interest)
  # interactome0 = integer pairs (interactome from which cluster0 comes)
  # this.interactome = integer pairs (noised interactome)
  
  # count edges in cluster of interest
  I0 = interactome0_int$protA %in% cluster0_int & interactome0_int$protB %in% cluster0_int
  N0 = sum(I)
  
  # count fraction of edges dropped 
  I = this.interactome$protA %in% cluster0_int & this.interactome$protB %in% cluster0_int
  tmp0 = interactome0_int$protA*10^6 + interactome0_int$protB
  tmp = this.interactome$protA*10^6 + this.interactome$protB
  ff = sum(tmp0[I0] %in% tmp[I]) / sum(I0)
  
  # count fraction of edges dropped for entire interactome
  FF = sum(tmp0 %in% tmp) / length(I0)
  
  E = ff / FF
  return(E)
}


geomacc = function(predComplex, refComplex){
  Na = length(refComplex)
  Nb = length(predComplex)
  
  Ni = numeric(Na)
  TT = matrix(nrow = Na, ncol = Nb)
  for (ii in 1:Na) {
    c1 = unlist(strsplit(refComplex[ii], ";"))
    Ni[ii] = length(c1)
    for (jj in 1:Nb) {
      c0 = unlist(strsplit(predComplex[jj], ";"))
      TT[ii,jj] = length(intersect(c1, c0));
    }
  }
  
  Trefmax = apply(TT,1,max)
  Sn = sum(Trefmax) / sum(Ni)
  
  Tpredmax = apply(TT,2,max)
  PPV = sum(Tpredmax) / sum(sum(TT))
  
  ga = sqrt(Sn * PPV)
  return(ga)
}


matchingratio = function(predComplex, refComplex){
  Na = length(refComplex)
  Nb = length(predComplex)
  
  TT = matrix(nrow = Na, ncol = Nb)
  for (ii in 1:Na) {
    c1 = unlist(strsplit(refComplex[ii], ";"))
    for (jj in 1:Nb) {
      c0 = unlist(strsplit(predComplex[jj], ";"))
      
      overlap = length(intersect(c1, c0)) ^2
      TT[ii,jj] = overlap / length(c1) / length(c0);
    }
  }
  
  if (Na<Nb) TT = t(TT)
  
  sortMatrix = matrix(nrow = nrow(TT), ncol = ncol(TT))
  for (ii in 1:ncol(TT)) {
    sortMatrix[,ii] = order(TT[,ii], decreasing = 1)
  }
  
  already_picked = numeric(nrow(TT))
  edges = numeric(ncol(TT))
  
  for (ii in 1:nrow(TT)) {
    tmp = numeric(ncol(TT))
    for (jj in 1:ncol(TT)) {
      tmp[jj] = TT[sortMatrix[ii,jj], jj]
    }
    Iorder = order(tmp, decreasing = 1)
    
    for (jj in 1:ncol(TT)) {
      Ipred = Iorder[jj]
      Iref = sortMatrix[ii, Ipred]
      
      if (!edges[Ipred]==0) next()
      
      if (already_picked[Iref]==0) {
        edges[Ipred] = Iref
        already_picked[Iref] = 1
      }
    }
  }
  
  mmr = 0;
  for (ii in 1:min(Na,Nb)) {
    if (edges[ii]>0) {
      mmr = mmr + TT[edges[ii], ii]
    }
  }
  mmr = mmr / Na
  
  return(mmr)
}



hiclust = function(ints, nclust){
  unqprots = unique(c(ints[,1], ints[,2]))
  nn = length(unqprots)
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints[,1], unqprots) # row
  I.col = match(ints[,2], unqprots) # column
  
  if (FALSE) { # long way
    # dummy dist object
    x = matrix(runif(nn * nn), nrow = nn, ncol=10)
    d.long = dist(x)
    attr(d.long, 'Upper') = T
    d.long[1:length(d.long)] = 1
    cc = 0
    for (ii in 1:nn) { # column
      ia = which(I.col==ii)
      for (jj in 1:nn) { # row
        if (ii>=jj) next
        cc = cc+1
        ib = which(I.row==jj)
        if (length(intersect(ia,ib))==1) {
          d.long[cc] = 0
        }
      }
    }
    hc = hclust(d.long)
    clusts = cutree(hc,5)
    
  } else { # short hard way
    I.fill = numeric(length(I.row))
    for (ii in 1:length(I.row)) {
      a = I.row[ii]
      b = I.col[ii]
      # ensure I.col < I.row, i.e. upper triangular
      if (I.row[ii] < I.col[ii]) {
        a = I.col[ii]
        b = I.row[ii]
      }
      
      I.fill[ii] = a - 1
      if (b>1) {
        colsum = 0
        for (jj in 1:(b-1)) {
          colsum = colsum + (nn-jj) - 1
        }
        I.fill[ii] = colsum + a - 1
      }
    }
    
    # dummy dist object
    x = matrix(runif(nn * 10), nrow = nn, ncol=10)
    d = dist(x)
    attr(d, 'Upper') = T
    d[1:length(d)] = 1
    d[I.fill] = 0
    hc = hclust(d)
    clusts = cutree(hc,nclust)
  }
  
  # compile `clusts` into lists of proteins
  unqclusts = unique(clusts)
  Nmembers = numeric(length(unqclusts))
  clusts.prots = character(length(unqclusts))
  for (ii in 1:length(unqclusts)) {
    I = clusts==unqclusts[ii]
    clusts.prots[ii] = paste(unqprots[I], collapse=";")
    Nmembers[ii] = sum(I)
  }
  clusts.prots = clusts.prots[Nmembers>=3]
  
  return(clusts.prots)
}


# hierachical
hierarch.edge.list.format = function(ints) {
  return(pam.edge.list.format(ints))
}

hierarch.cluster.format = function(tmp, unqprots) {
  clusts = character()
  unqclusts = unique(tmp)
  for (ii in 1:length(unqclusts)) {
    I = tmp == unqclusts[ii]
    clusts[ii] = paste(unqprots[I], collapse = ";")
  }
  clusts = clusts[!clusts==""]
  clusts = clusts[!is.na(clusts)]
  return(clusts)
}


pamclust = function(ints, nclust){
  unqprots = unique(c(ints[,1], ints[,2]))
  nn = length(unqprots)
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints[,1], unqprots) # row
  I.col = match(ints[,2], unqprots) # column
  
  unqprots = unique(c(ints[,1], ints[,2]))
  nn = length(unqprots)
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints[,1], unqprots) # row
  I.col = match(ints[,2], unqprots) # column
  
  I.fill = numeric(length(I.row))
  for (ii in 1:length(I.row)) {
    a = I.row[ii]
    b = I.col[ii]
    # ensure I.col < I.row, i.e. upper triangular
    if (I.row[ii] < I.col[ii]) {
      a = I.col[ii]
      b = I.row[ii]
    }
    
    I.fill[ii] = a - 1
    if (b>1) {
      colsum = 0
      for (jj in 1:(b-1)) {
        colsum = colsum + (nn-jj) - 1
      }
      I.fill[ii] = colsum + a - 1
    }
  }
  
  # dummy dist object
  x = matrix(runif(nn * 10), nrow = nn, ncol=10)
  d = dist(x)
  attr(d, 'Upper') = T
  d[1:length(d)] = 1
  d[I.fill] = 0
  if (ncol(ints)==3) d[I.fill] = ints[,3]
  
  clusts = pam(d, nclust)
  
  # distmat = as.matrix(matrix(nrow=nn, ncol=nn))
  # distmat[is.na(distmat)] = 1
  # for (ii in 1:length(I.row)) {
  #   distmat[I.row[ii], I.col[ii]] = runif(1)*.0001
  #   distmat[I.col[ii], I.row[ii]] = runif(1)*.0001
  # }
  
  # MDS
  #fit = cmdscale(distmat,eig=TRUE, k=5)
  #fit = isoMDS(distmat, k=2)
  
  # PAM
  #clusts = pam(distmat, nclust)
  
  # fastkmed
  #clusts = fastkmed(distmat, ncluster = nclust, iterate = 50)
  
  # compile `clusts` into lists of proteins
  unqclusts = unique(clusts$clustering)
  Nmembers = numeric(length(unqclusts))
  clusts.prots = character(length(unqclusts))
  for (ii in 1:length(unqclusts)) {
    I = clusts$cluster==unqclusts[ii]
    clusts.prots[ii] = paste(unqprots[I], collapse=";")
    Nmembers[ii] = sum(I)
  }
  clusts.prots = clusts.prots[Nmembers>=3]
  
  return(clusts.prots)
}



shufflecorum = function(ints.corum, ff){
  unqprots = sort(unique(c(ints.corum[,1], ints.corum[,2])))
  unqprots = unqprots[!unqprots==""]
  
  #
  ints.shuffle = ints.corum
  N.replace = round(nrow(ints.corum) * ff)
  I.replace = sample(nrow(ints.corum), N.replace)
  
  # indices of unshuffled
  ia0 = match(ints.corum[,1], unqprots)
  ib0 = match(ints.corum[,2], unqprots)

  # indices of shuffled
  ia = ia0
  ib = ib0
  ia[I.replace] = sample(length(unqprots), N.replace, replace = T)
  ib[I.replace] = sample(length(unqprots), N.replace, replace = T)
  
  # ensure no self-interactions
  while (sum(ia==ib)>0) {
    I = ia==ib
    ib[I] = sample(length(unqprots), sum(I), replace = T)
  }
  
  # ensure no duplicates
  while (sum(duplicated(c(paste(ia,ib))))>0) {
    I = duplicated(c(paste(ia,ib)))
    ib[I] = sample(length(unqprots), sum(I), replace = T)
  }
  
  #
  ints.shuffle[,1] = unqprots[ia]
  ints.shuffle[,2] = unqprots[ib]
  if (ncol(ints.corum)==3) {
    score.shuffle = ints.corum[,3]
    score.shuffle[I.replace] = sample(scores, N.replace)
    ints.shuffle[,3] = score.shuffle
  }
  
  # quality control: make sure you shuffled N.replace interactions
  N.diff = sum(!ia0==ia | !ib0==ib)
  if (abs(N.replace - N.diff)>5) print(paste("shuffling missed", N.replace-N.diff, "interactions"))
  
  # # sort protein pairs alphabetically
  # for (ii in 1:length(I)) {
  #   tmp = sort(c(ia[ii], ib[ii]))
  #   ia[ii] = tmp[1]
  #   ib[ii] = tmp[2]
  # }
  
  return(ints.shuffle) 
}

addcorum = function(ints.corum, ff){
  unqprots = sort(unique(c(ints.corum[,1], ints.corum[,2])))
  unqprots = unqprots[!unqprots==""]
  
  N.add = round(nrow(ints.corum) * ff)
  
  # N random pairs
  protA = sample(unqprots, N.add, replace = T)
  protB = sample(unqprots, N.add, replace = T)
  
  I.bad = 1
  while (sum(I.bad)>0) {
    # sort A<B
    ia = match(protA, unqprots)
    ib = match(protB, unqprots)
    I = ia>ib
    protA0 = protA
    protB0 = protB
    protA[I] = protB0[I]
    protB[I] = protA0[I]
    
    # replace any A==B or AB%in%corum
    I.bad = protA==protB | (paste(protA,protB) %in% paste(ints.corum[,1], ints.corum[,2]))
    protA[I.bad] = sample(unqprots, sum(I.bad), replace = T)
    protB[I.bad] = sample(unqprots, sum(I.bad), replace = T)
  }
  
  ints.add = data.frame(protA = protA, protB = protB, stringsAsFactors = F)
  if (ncol(ints.corum)==2) {
    names(ints.add) = names(ints.corum)
  } else if (ncol(ints.corum)==3) {
    ints.add$score = sample(ints.corum[,3], nrow(ints.add))
    names(ints.add) = names(ints.corum)
  }
  ints.shuffle = rbind(ints.corum, ints.add)
  return(ints.shuffle)
}

removecorum = function(ints.corum, ff){
  N.keep = round(nrow(ints.corum) * (1 - ff))
  I.keep = sample(nrow(ints.corum), N.keep)
  ints.corum = ints.corum[I.keep, ]
  return(ints.corum)
}

consensus.adjmat = function(this.cluster, clusters){
  unqIters = unique(clusters$iter)
  
  this.cluster = unlist(strsplit(this.cluster, ";"))
  this.size = length(this.cluster)
  this.size = formatC(this.size, width=3, flag="0")
  if (length(this.cluster)>150) next
  
  # find its best match in the other iters
  tmp = numeric(10)
  bestMatches = character(length(unqIters))
  for (mm in 1:length(unqIters)) {
    set0 = clusters$cluster[clusters$iter==unqIters[mm]]
    JJ = numeric(length(set0))
    for (uu in 1:length(set0)) {
      that.cluster = unlist(strsplit(set0[uu], ";"))
      JJ[uu] = length(intersect(this.cluster, that.cluster)) / 
        length(unique(c(this.cluster, that.cluster)))
    }
    tmp[mm] = max(JJ)
    I.match = which.max(JJ)
    bestMatches[mm] = set0[I.match]
  }
  
  # make consensus adjacency matrix
  allProts = unique(c(this.cluster, unlist(lapply(bestMatches, strsplit, ";"))))
  # put this.cluster in the middle
  tmp = allProts[!allProts %in% this.cluster]
  if (length(tmp)==1) {
    allprots = c(this.cluster, tmp)
  } else if (length(tmp)>1) {
    nn = round(length(tmp)/2)
    allProts = c(tmp[1:nn], this.cluster, tmp[(nn+1):length(tmp)])
  }
  adjmat = matrix(numeric(length(allProts)^2), nrow=length(allProts), ncol=length(allProts))
  df.adjmat = data.frame(prots = character(0),
                         variable = character(0),
                         value = numeric(0), iter=numeric(0), stringsAsFactors = F)
  bestMatches = c(paste(this.cluster,collapse=";"), bestMatches)
  Ar = numeric(length(bestMatches))
  for (mm in 1:length(bestMatches)) {
    that.cluster = unlist(strsplit(bestMatches[mm], ";"))
    for (uu in 1:length(that.cluster)) {
      ia = which(allProts %in% that.cluster[uu])
      for (vv in 1:length(that.cluster)) {
        if (uu>=vv) next
        ib = which(allProts %in% that.cluster[vv])
        adjmat[ia,ib] = adjmat[ia,ib]+1
        adjmat[ib,ia] = adjmat[ib,ia]+1
      }
    }
    
    Ar[mm] = calcA(this.cluster, paste(that.cluster, collapse=";"))
    df.tmp = as.data.frame(adjmat)
    df.tmp$prots = allProts
    df.tmp2 = melt(df.tmp, id.vars="prots")
    df.tmp2$iter = mm
    df.tmp2$value = df.tmp2$value * length(bestMatches) / mm
    df.adjmat = rbind(df.adjmat, df.tmp2)
  }
  df.adjmat$prots = factor(df.adjmat$prots, levels = allProts)
  df.adjmat$variable = levels(df.adjmat$prots)[match(df.adjmat$variable, levels(df.adjmat$variable))]
  df.adjmat$variable = factor(df.adjmat$variable, levels = allProts)
  
  return(df.adjmat)
}


calcMIz = function(X,Y, iterMax, error.bars=F) {
  
  nn = length(X)
  zz = numeric(iterMax)
  for (iter in 1:iterMax) {
    # shuffle one of the labels
    I = sample(nn, nn)
    zz[iter] = mutinformation(X[I],Y)
  }
  mi = mutinformation(X,Y)
  miz = (mi - mean(zz)) / sd(zz)
  
  if (error.bars) {
    miz.error = numeric(100)
    for (ii in 1:100) {
      this.zz = sample(zz, round(iterMax/2))
      miz.error[ii] = (mi - mean(this.zz)) / sd(this.zz)
    }
    miz = (mi - mean(this.zz)) / sd(zz)
    return(c(quantile(miz.error, .025), mean(miz.error), quantile(miz.error, .975)))
  } else {
    return(miz)
  }
  
}


# make your own melt function
meltmat = function(mat, id.vars) {
  nmes = names(mat)
  I = !nmes %in% id.vars
  mat.to.melt = mat[,I]
  
  nn = nrow(mat.to.melt)
  nn2 = nn^2
  
  df = data.frame(x=numeric(nn2), y=numeric(nn2), value=numeric(nn2), proteins=numeric(nn2))
  cc = 0
  for (ii in 1:nrow(mat.to.melt)) {
    for (jj in 1:ncol(mat.to.melt)) {
      cc = cc+1
      df$x[cc] = ii
      df$y[cc] = jj
      df$value[cc] = mat.to.melt[ii,jj]
      df$proteins[cc] = mat$proteins[ii]
    }
  }
  
  return(df)
}


binarize.corum = function(corum) {
  
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
    ints.corum[,1][I] = pairs.prots[,1]
    ints.corum[,2][I] = pairs.prots[,2]
    cc = cc+length(I)
  }
  ints.corum = ints.corum[1:cc,]
  ints.corum = distinct(ints.corum)
  
  return(ints.corum)
}


netdiff = function(g1, g2) {
  
  df = data.frame(edge.gain = numeric(1), 
                  edge.loss = numeric(1), 
                  node.gain = numeric(1),
                  node.loss = numeric(1),
                  En1n1 = numeric(1),
                  En1n2 = numeric(1),
                  En2n2 = numeric(1))
  
  prots1 = names(V(g1))
  prots2 = names(V(g2))
  ints1 = attr(E(g1), "vnames")
  ints2 = attr(E(g2), "vnames")
  
  df$edge.gain = sum(! ints2 %in% ints1)
  df$edge.loss = sum(! ints1 %in% ints2)
  df$node.gain = sum(! prots2 %in% prots1)
  df$node.loss = sum(! prots1 %in% prots2)
  
  n2 = prots2[!prots2 %in% prots1] # nodes only in net2
  
  tmp = sapply(ints2, FUN = function(x) strsplit(x, "|", fixed=T))
  na = unlist(sapply(tmp, "[", 1))
  nb = unlist(sapply(tmp, "[", 2))
  i.new = !ints2 %in% ints1
  
  df$En1n1 = sum(!na[i.new]%in%n2 & !nb[i.new]%in%n2) # edge between original nodes
  df$En2n2 = sum(na[i.new]%in%n2 & nb[i.new]%in%n2) # edge between new nodes
  df$En1n2 = sum(i.new) - df$En1n1 - df$En2n2
  
  return(df)
}



# dependencies
require(igraph)
require(readr)
require(Matrix)
require(cluster)

# functions

calcJ = function(this.cluster, clusters) {
  # J, Ji
  # assignment reproducibility
  # essentially the maximum Jaccard index
  #
  # this.cluster = single string, with semicolon-delimited IDs
  # clusters = vector of strings, all semicolon-delimited IDs
  
  this.cluster = unlist(strsplit(this.cluster, ";"))
  JJ = numeric(length(clusters))
  for (ii in 1:length(clusters)) {
    that.cluster = unlist(strsplit(clusters[ii], ";"))
    JJ[ii] = length(intersect(this.cluster, that.cluster)) / 
      length(unique(c(this.cluster, that.cluster)))
  }
  Ji = max(JJ, na.rm=T)
  tmp = list(Ji = Ji, best.cluster = clusters[which.max(JJ)])
  
  return(tmp)
}

shufflecorum = function(ints.corum, ff){
  unqprots = sort(unique(c(ints.corum[,1], ints.corum[,2])))
  unqprots = unqprots[!unqprots==""]
  
  #
  ints.shuffle = ints.corum
  N.replace = round(nrow(ints.corum) * ff)
  I.replace = sample(nrow(ints.corum), N.replace)
  
  # indices of unshuffled
  ia0 = match(ints.corum[,1], unqprots)
  ib0 = match(ints.corum[,2], unqprots)
  
  # indices of shuffled
  ia = ia0
  ib = ib0
  ia[I.replace] = sample(length(unqprots), N.replace, replace = T)
  ib[I.replace] = sample(length(unqprots), N.replace, replace = T)
  
  # ensure no self-interactions
  while (sum(ia==ib)>0) {
    I = ia==ib
    ib[I] = sample(length(unqprots), sum(I), replace = T)
  }
  
  #
  ints.shuffle[,1] = unqprots[ia]
  ints.shuffle[,2] = unqprots[ib]
  
  # quality control: make sure you shuffled N.replace interactions
  N.diff = sum(!ia0==ia | !ib0==ib)
  #if (abs(N.replace - N.diff)>5) print(paste("shuffling missed", N.replace-N.diff, "interactions"))
  
  return(ints.shuffle) 
}


pam.edge.list.format = function(ints) {
  unqprots = unique(c(ints[,1], ints[,2]))
  nn = length(unqprots)
  scores = rep(0, nrow(ints))
  if (ncol(ints)==3) scores = ints[,3]
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints[,1], unqprots) # row
  I.col = match(ints[,2], unqprots) # column
  
  unqprots = unique(c(ints[,1], ints[,2]))
  nn = length(unqprots)
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints[,1], unqprots) # row
  I.col = match(ints[,2], unqprots) # column
  
  I.fill = numeric(length(I.row))
  for (ii in 1:length(I.row)) {
    a = I.row[ii]
    b = I.col[ii]
    # ensure I.col < I.row, i.e. upper triangular
    if (I.row[ii] < I.col[ii]) {
      a = I.col[ii]
      b = I.row[ii]
    }
    
    I.fill[ii] = a - 1
    if (b>1) {
      colsum = 0
      for (jj in 1:(b-1)) {
        colsum = colsum + (nn-jj) - 1
      }
      I.fill[ii] = colsum + a - 1
    }
  }
  
  # dummy dist object
  x = matrix(stats::runif(nn * 10), nrow = nn, ncol=10)
  d = stats::dist(x)
  attr(d, 'Upper') = T
  d[1:length(d)] = 1
  d[I.fill] = scores
  
  return(d)
}

pam.cluster.format = function(clusts, unqprots) {
  # compile `clusts` into lists of proteins
  unqclusts = unique(clusts$clustering)
  Nmembers = numeric(length(unqclusts))
  clusts.prots = character(length(unqclusts))
  for (ii in 1:length(unqclusts)) {
    I = clusts$cluster==unqclusts[ii]
    clusts.prots[ii] = paste(unqprots[I], collapse=";")
    Nmembers[ii] = sum(I)
  }
  #clusts.prots = clusts.prots[Nmembers>=3]
  #clusts.prots = as.list(clusts.prots)
  
  return(clusts.prots)
}


mcl.edge.list.format = function(ints.corum) {
  G = igraph::graph.data.frame(ints.corum,directed=FALSE)
  A = igraph::as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)
}

mcl.cluster.format = function(tmp, unqnodes) {
  tmp = tmp$Cluster
  clusts = character()
  unqclusts = unique(tmp)
  for (ii in 1:length(unqclusts)) {
    I = tmp == unqclusts[ii]
    #if (sum(I)<3) next
    clusts[ii] = paste(unqnodes[I], collapse = ";")
  }
  clusts = clusts[!clusts==""]
  clusts = clusts[!is.na(clusts)]
  return(clusts)
}

# hierachical
hierarch.edge.list.format = function(ints) {
  return(pam.edge.list.format(ints))
}

hierarch.cluster.format = function(tmp, unqprots) {
  clusts = character()
  unqclusts = unique(tmp)
  for (ii in 1:length(unqclusts)) {
    I = tmp == unqclusts[ii]
    clusts[ii] = paste(unqprots[I], collapse = ";")
  }
  clusts = clusts[!clusts==""]
  clusts = clusts[!is.na(clusts)]
  return(clusts)
}

get.full.analysis = function(lazyload = T) {
  # get full grid
  #   (mcl, co_mcl, co, pam, walktrap, hierarch, mcode, louvain, leiden) x (corum, drugbank, emailEU, biogrid, collins2007)
  
  if (lazyload==F) {
    # (mcl, co_mcl, co, pam, walktrap) x (corum, drugbank, emailEU)
    sf = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/data/fig2B_v06.Rda"
    load(sf)
    
    # (mcode, louvain, leiden, hierarch,)x(corum, drugbank, emailEU) + 
    #     (mcl, co_mcl, co, pam, walktrap, hierarch, mcode, louvain, leiden)x(biogrid, collins2007)
    fns = list.files("../data/data 2/clusters/", pattern = "*.txt", full.names = T)
    bad = rep(NA, 1000)
    for (ii in 1:length(fns)) {
      tmp = as.data.frame(read_tsv(fns[ii]))
      if (nrow(tmp)==0) bad[ii] = fns[ii]
      if (!length(names(tmp))==4) next
      names(tmp) = c("cluster", "algorithm", "noise_mag", "network")
      tmp$Ji1 = rep(NA, nrow(tmp))
      tmp$iter = rep(NA, nrow(tmp))
      tmp$cluster.size = unlist(lapply(sapply(tmp$cluster, strsplit, ";"), length))
      tmp = tmp[names(Ji)]
      Ji = rbind(Ji, tmp)
    }
    
    # clean up
    Ji$algorithm[Ji$algorithm=="CO"] = "co"
    Ji$algorithm[Ji$algorithm=="MCL"] = "mcl"
    Ji$network[Ji$network=="ChCh-Miner_durgbank-chem-chem.tsv"] = "DrugBank"
    Ji$network[Ji$network=="BioPlex_293T_Network_10K_Dec_2019.tsv"] = "BioPlex"
    Ji$network[Ji$network=="email-Eu-core.txt"] = "email-Eu"
    Ji$network[Ji$network=="corum_pairwise.txt"] = "CORUM"
    Ji$network[Ji$network=="BIOGRID-ALL-3.5.186.tab3.txt"] = "BIOGRID"
    Ji$network[Ji$network=="PE_and_conf_scores_01-11-08.txt"] = "Collins2007"
    
    unqnet = unique(Ji$network)
    unqalg = unique(Ji$algorithm)
    unqnoise = sort(unique(Ji$noise))
    for (ii in 1:length(unqnet)) {
      for (jj in 1:length(unqalg)) {
        I0 = which(Ji$network==unqnet[ii] & Ji$algorithm==unqalg[jj] & Ji$noise_mag==0)
        for (kk in 1:length(unqnoise)) {
          if (kk==1) {
            Ji$Ji1[I0] = 1
            next
          }
          
          I1 = which(Ji$network==unqnet[ii] & Ji$algorithm==unqalg[jj] & Ji$noise_mag==unqnoise[kk])
          if (sum(I1)==0) print(paste("no data for", unqnet[ii], unqalg[jj], unqnoise[kk]))
          next
          cluster0 = Ji$cluster[I0]
          cluster = Ji$cluster[I1]
          Jii = numeric(length(I1))
          for (mm in 1:length(cluster)) {
            Ji[mm] = calcA(cluster[mm], cluster0)
          }
          Ji$Ji1[I1] = Ji
        }
      }
    }
  } else {
    fn = ""
    load(fn)
  }
  
  return(Ji)
}

get.addremove.analysis = function(lazyload = T) {
  
  if (lazyload==F) {
    fns = list.files("../data/clusters_add_remove//", pattern = "*.txt", full.names = T)
    Ji3 = as.data.frame(read_tsv(fns[1]))
    for (ii in 2:length(fns)) {
      tmp = as.data.frame(read_tsv(fns[ii]))
      if (nrow(tmp)==0) bad[ii] = fns[ii]
      Ji3 = rbind(Ji3, tmp)
    }

    # calculate Ji
    fn.save = "../data/mcp_mcode_addremove.txt"
    Ji3$Ji = rep(NA, nrow(Ji3))
    unqalg = unique(Ji3$algorithm)
    unqadd = unique(Ji3$add_mag)
    unqremove = unique(Ji3$remove_mag)
    for (ii in 1:length(unqalg)) {
      ia = Ji3$algorithm == unqalg[ii]
      cluster0 = Ji3$cluster[ia & Ji3$add_mag==0 & Ji3$remove_mag==0]
      for (jj in 1:length(unqadd)) {
        ib = Ji3$add_mag==unqadd[jj]
        for (kk in 1:length(unqremove)) {
          ic = Ji3$remove_mag==unqremove[kk]
          II = ia & ib & ic
          cluster = Ji3$cluster[II]
          
          print(paste(unqalg[ii], " add-", unqadd[jj], " remove-", unqremove[kk], sep=""))
          
          Jii = numeric(sum(II))
          for (mm in 1:length(cluster)) {
            Jii[mm] = calcA(cluster[mm], cluster0)
          }
          Ji3$Ji[II] = Jii
          
          write_tsv(Ji3, path = fn.save)
        }
      }
    }
    
  } else {
    fn = "../data/mcp_addremove.txt"
    Ji3 = as.data.frame(read_tsv(fn))
    Ji3 = Ji3[!Ji3$algorithm == "mcode", ]
    
    fn = "../data/mcp_mcode_addremove.txt"
    tmp =  as.data.frame(read_tsv(fn))
    Ji3 = rbind(Ji3, tmp)
  }
  
  return(Ji3)
}


