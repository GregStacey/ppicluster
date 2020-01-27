
# dependencies


# functions

calcA2 = function(this.cluster, clusters) {
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
  ia0 = match(ints.corum$protA, unqprots)
  ib0 = match(ints.corum$protB, unqprots)
  
  # defined original pairs, possible new pairs
  tmp = t(combn(1:length(unqprots), 2))
  all = paste(tmp[,1], tmp[,2])
  orig = paste(ia0, ib0)
  poss = all[!all %in% orig]
  
  # indices of shuffled
  ia = ia0
  ib = ib0
  I.new = sample(poss, length(I.replace))
  ia[I.replace] = as.numeric(sapply(sapply(I.new, strsplit, " "), "[", 1))
  ib[I.replace] = as.numeric(sapply(sapply(I.new, strsplit, " "), "[", 2))
  #ia[I.replace] = sample(length(unqprots), N.replace, replace = T)
  #ib[I.replace] = sample(length(unqprots), N.replace, replace = T)
  
  # ensure no self-interactions
  while (sum(ia==ib)>0) {
    I = ia==ib
    ib[I] = sample(length(unqprots), sum(I), replace = T)
  }
  
  # ensure no duplicates
  while (sum(duplicated(c(paste(ia,ib))))>0) {
    I = which(duplicated(c(paste(ia,ib))))
    for (ii in 1:length(J)) {
      good = (1:length(unqprots))
      good = good[!good %in% which(ia0[ib[I[ii]]]==1)]
      ib[I[ii]] = sample(length(unqprots), sum(I), replace = T)
    }
  }
  
  #
  ints.shuffle$protA = unqprots[ia]
  ints.shuffle$protB = unqprots[ib]
  
  # quality control: make sure you shuffled N.replace interactions
  N.diff = sum(!ia0==ia | !ib0==ib)
  if (abs(N.replace - N.diff)>5) print(paste("shuffling missed", N.replace-N.diff, "interactions"))
  
  return(ints.shuffle) 
}

calc.reproducibility = function(this.cluster, clusters){
  unqIters = unique(clusters$iter)
  
  this.cluster = unlist(strsplit(this.cluster, ";"))
  this.size = length(this.cluster)
  this.size = formatC(this.size, width=3, flag="0")
  
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
  nn = round(length(tmp)/2)
  allProts = c(tmp[1:nn], this.cluster, tmp[(nn+1):length(tmp)])
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
  
  return(df.adjmat)
}


pam.edge.list.format = function(ints) {
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
  clusts.prots = clusts.prots[Nmembers>=3]
  clusts.prots = as.list(clusts.prots)
  
  return(clusts.prots)
}

# mcl
# requires passing cluster.format the unqprots from the edge list
# this is so clusters, e.g. "1;2;3" are matched to proteins
mcl.edge.list.format = function(ints.corum) {
  G = graph.data.frame(ints.corum,directed=FALSE)
  A = as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)
}

mcl.cluster.format = function(tmp, unqprots) {
  clusts = character()
  unqclusts = unique(tmp)
  for (ii in 1:length(unqclusts)) {
    clusts[ii] = paste(unqprots[tmp == unqclusts[ii]], collapse = ";")
  }
  return(clusts)
}
