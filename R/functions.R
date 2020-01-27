
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
require(VennDiagram)
require(infotheo)
require(igraph)
require(dils)
#require(MCL)
require(hbm)
mcl = hbm::mcl
require(ontologyIndex)
require(flavin)
source("clusterone_java.R")


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


hairball = function(adjmat, this.cluster, allProts){
  # convert to graph
  g = graph_from_adjacency_matrix(adjmat, mode = 'undirected', weighted=T)
  I.orig = (allProts %in% this.cluster)
  
  # add information to each node 
  nodes = allProts
  palette = colorRamp(ggsci::pal_gsea()(n = 12))(seq(100) / 100)
  colors = palette[E(g)$weight * 255]
  V(g)$color[I.orig] = "black"
  V(g)$color[!I.orig] = "white"
  V(g)$size[I.orig] = .25
  V(g)$size[!I.orig] = 0
  E(g)$alpha = E(g)$weight
  E(g)$size = E(g)$weight
  I.orig = which(I.orig)
  I.orig.pairs = numeric(length(I.orig) * (length(I.orig)-1) / 2)
  cc = 0
  for (ii in 1:length(I.orig)) {
    for (jj in 1:length(I.orig)) {
      if (ii>=jj) next
      cc = cc+1
      I.orig.pairs[cc] = I.orig[ii]
      cc = cc+1
      I.orig.pairs[cc] = I.orig[jj]
    }
  }
  E(g)$color = 0
  E(g, P = I.orig.pairs)$color = 1
  
  layout = layout.lgl(g)
  N = ggnetwork(g, layout = layout)
  N$cluster0 = as.numeric(N$vertex.names %in% I.orig)
  # add two extra line to prevent normalizing
  N$size.x = N$size.x / 4
  x = N[nrow(N),]
  N = rbind(N,x)
  N = rbind(N,x)
  iend = nrow(N)
  N[iend,c("weight","size.x","color.y","alpha","cluster0")] = 1
  N[iend-1,c("weight","size.x","color.y","alpha","cluster0")] = 0
  N[iend,c("x","y","xend","yend")] = 0
  N$x[iend] = 0.125
  
  p = ggplot(N, aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_edges(aes(size=size.x, color=color.y, alpha=alpha)) + 
    geom_nodes(aes(color=cluster0, size = cluster0)) +
    theme_void() + theme(legend.position="none")
  return(p)
}


hiclust = function(ints, nclust){
  unqprots = unique(c(ints$protA, ints$protB))
  nn = length(unqprots)
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints$protA, unqprots) # row
  I.col = match(ints$protB, unqprots) # column
  
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




pamclust = function(ints, nclust){
  unqprots = unique(c(ints$protA, ints$protB))
  nn = length(unqprots)
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints$protA, unqprots) # row
  I.col = match(ints$protB, unqprots) # column

  unqprots = unique(c(ints$protA, ints$protB))
  nn = length(unqprots)
  
  # create adjacency matrix in the style of `dist` object
  # MAKE SURE YOU'RE GETTING THIS RIGHT!!!
  I.row = match(ints$protA, unqprots) # row
  I.col = match(ints$protB, unqprots) # column
  
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
  ia0 = match(ints.corum$protA, unqprots)
  ib0 = match(ints.corum$protB, unqprots)
  
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
  ints.shuffle$protA = unqprots[ia]
  ints.shuffle$protB = unqprots[ib]
  
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
    I.bad = protA==protB | (paste(protA,protB) %in% paste(ints.corum$protA, ints.corum$protB))
    protA[I.bad] = sample(unqprots, sum(I.bad), replace = T)
    protB[I.bad] = sample(unqprots, sum(I.bad), replace = T)
  }
  
  ints.add = data.frame(protA = protA, protB = protB, stringsAsFactors = F)
  ints.corum = rbind(ints.corum, ints.add)
  return(ints.corum)
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
    ints.corum$protA[I] = pairs.prots[,1]
    ints.corum$protB[I] = pairs.prots[,2]
    cc = cc+length(I)
  }
  ints.corum = ints.corum[1:cc,]
  ints.corum = distinct(ints.corum)
  
  return(ints.corum)
}
