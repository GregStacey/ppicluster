
# TO DO
# 
# 1. Include permutation-based significance for repJ numbers
#    This will be very time intensive, so make it optional.



#' Perturb clusters
#' 
#' Network to be clustered is formated as a data frame edge list (edge.list). Users can 
#' pass a custom function edge.list.format to handle different formats required by
#' clustering.algorithm. edge.list.format must transform a data frame edge list (two 
#' columns unweighted, or three columsn weighted) into the required format.
#' 
#' Assumes clustering.algorithm returns a list of clustered nodes. If clustering.algorithm
#' returns a different format, users can pass a custom function cluster.format that
#' transforms the output of clustering.algorithm into a list of clustered nodes.
#' 
#' @param edge.list data frame with two or three columns, with first two columns providing 
#' the list of edges and an optional third column giving edge weight.
#' @param clustering.algorithm function responsible for clustering.
#' @param noise vector with values between 0 and 1 specifying the amount(s) of noise to 
#' add to the network.
#' @param iters positive integer specifying number of iterations.
#' @param edge.list.format optional function that transforms edge list into format required
#' by clustering.algorithm.
#' @param cluster.format optional function that transforms output returned by 
#' clustering.algorithm into a data frame with three or four columns, where columns one 
#' and two are an edge list, column three is an id specifying the cluster assignment, and 
#' @return data frame containing all clusters and their density score.
#' @examples
#' first example
#' 
#' second example

clust.perturb = function(network, 
                         clustering.algorithm, 
                         noise = 0.5, 
                         iters = 100, 
                         edge.list.format = NULL,
                         cluster.format = NULL) {
  
  #print(clustering.algorithm)
  #print(edge.list.format)
  #print(cluster.format)
  
  # cluster without noise
  cc = 0
  network.input = network
  if (!is.null(edge.list.format)) network.input = edge.list.format(network)
  tmp = clustering.algorithm(network.input)
  if (!is.null(cluster.format)) tmp = cluster.format(tmp)
  # store clusters
  clusters0 = data.frame(cluster = character(length(tmp)),
                         reproducibility.J = rep(NA, length(tmp)),
                         all.iterations.repJ = rep(NA, length(tmp)),
                         #best.matches = rep(NA, length(tmp)),
                         stringsAsFactors = F)
  for (jj in 1:length(tmp)) clusters0$cluster[jj] = paste(tmp[[jj]], collapse=";")
  rm(tmp)
  str(clusters0)
  
  # cluster with noise
  clusters.noise = data.frame(iter = numeric(10^6),
                              noise_mag = numeric(10^6),
                              cluster = character(10^6), stringsAsFactors = F)
  cc = 0
  for (iter in 1:iters) {
    for (ii in 1:length(noise)) {
      print(paste("clustering iter",iter,"at noise=", noise[ii]))
      
      # add noise to network
      ints.shuffle = shufflecorum(network, noise[ii])
      
      # transform network to required format (if needed)
      if (!is.null(edge.list.format)) ints.shuffle = edge.list.format(ints.shuffle)
      
      # cluster
      these.clusters = clustering.algorithm(ints.shuffle)
      
      # transform clusters to list (if needed)
      if (!is.null(cluster.format)) these.clusters = cluster.format(these.clusters)
      
      # store these.clusters
      for (jj in 1:length(these.clusters)) {
        cc = cc+1
        clusters.noise$iter[cc] = iter
        clusters.noise$noise_mag[cc] = noise[ii]
        clusters.noise$cluster[cc] = paste(these.clusters[[jj]], collapse=";")
      }
    }
  }
  clusters.noise = clusters.noise[1:cc,]
  
  
  # calculate J for every cluster in clusters0
  for (ii in 1:nrow(clusters0)) {
    print(paste("cluster", ii))
    tmp.j = numeric(iters * length(noise))
    best.match = character(iters * length(noise))
    cc = 0
    for (jj in 1:iters) {
      for (kk in 1:length(noise)) {
        noise.clusters = clusters.noise$cluster[clusters.noise$iter==jj & 
                                                  clusters.noise$noise_mag==noise[kk]]
        cc = cc+1
        tmp = calcA(clusters0$cluster[ii], noise.clusters)
        tmp.j[cc] = tmp[[1]]
        #best.match[ii] = tmp[[2]]
      }
    }
    
    clusters0$reproducibility.J[ii] = mean(tmp.j, na.rm=T)
    clusters0$all.iterations.repJ[ii] = paste(tmp.j, collapse = ";")
    #clusters0$best.matches[ii] = paste(best.match, collapse = " | ")
  }
  
  return(clusters0)
}



# the same as clust.perturb, except
#   - also returns the number of shuffled edges per cluster0 in each iter
#   - also returns the best matches in each iter
#   - also returns best50, best 75
clust.perturb2 = function(network, 
                          clustering.algorithm, 
                          noise = 0.5, 
                          iters = 100, 
                          edge.list.format = NULL,
                          cluster.format = NULL) {
  
  # cluster without noise
  cc = 0
  network.input = network
  if (!is.null(edge.list.format)) network.input = edge.list.format(network)
  tmp = clustering.algorithm(network.input)
  if (!is.null(cluster.format)) tmp = cluster.format(tmp)
  # store clusters
  clusters0 = data.frame(cluster = character(length(tmp)),
                         n.removed = rep(NA, length(tmp)),
                         reproducibility.J = rep(NA, length(tmp)),
                         all.iterations.repJ = rep(NA, length(tmp)),
                         best.match = rep(NA, length(tmp)),
                         best50 = rep(NA, length(tmp)),
                         best75 = rep(NA, length(tmp)),
                         stringsAsFactors = F)
  for (jj in 1:length(tmp)) clusters0$cluster[jj] = paste(tmp[[jj]], collapse=";")
  rm(tmp)
  str(clusters0)
  
  # cluster with noise
  n.changed = matrix(nrow = iters*length(noise), ncol = nrow(clusters0))
  clusters.noise = data.frame(iter = numeric(10^6),
                              noise_mag = numeric(10^6),
                              cluster = character(10^6), stringsAsFactors = F)
  cc = 0
  cc2 = 0
  for (iter in 1:iters) {
    for (ii in 1:length(noise)) {
      print(paste("clustering iter",iter,"at noise=", noise[ii]))
      cc2 = cc2+1
      
      # add noise to network
      ints.shuffle = shufflecorum(network, noise[ii])
      i.removed = !paste(network$protA,network$protB,sep="-") %in% paste(ints.shuffle$protA,ints.shuffle$protB,sep="-")
      ints.removed = network[i.removed,]
      ia = sapply(ints.removed$protA, FUN = function(x) grepl(x, clusts))
      ib = sapply(ints.removed$protB, FUN = function(x) grepl(x, clusts))
      n.changed[cc2,] = rowSums(ia & ib, na.rm=T)
      
      # connect ints.shuffle to clusters0
      
      # transform network to required format (if needed)
      if (!is.null(edge.list.format)) ints.shuffle = edge.list.format(ints.shuffle)
      
      # cluster
      these.clusters = clustering.algorithm(ints.shuffle)
      
      # transform clusters to list (if needed)
      if (!is.null(cluster.format)) these.clusters = cluster.format(these.clusters)
      
      # store these.clusters
      for (jj in 1:length(these.clusters)) {
        cc = cc+1
        clusters.noise$iter[cc] = iter
        clusters.noise$noise_mag[cc] = noise[ii]
        clusters.noise$cluster[cc] = paste(these.clusters[[jj]], collapse=";")
      }
    }
  }
  clusters.noise = clusters.noise[1:cc,]
  
  
  # calculate J for every cluster in clusters0
  clusters0$best50 = character(nrow(clusters0))
  for (ii in 1:nrow(clusters0)) {
    print(paste("cluster", ii))
    tmp.j = numeric(iters * length(noise))
    best.match = character(iters* length(noise))
    cc = 0
    for (jj in 1:iters) {
      for (kk in 1:length(noise)) {
        noise.clusters = clusters.noise$cluster[clusters.noise$iter==jj & 
                                                  clusters.noise$noise_mag==noise[kk]]
        cc = cc+1
        tmp = calcA2(clusters0$cluster[ii], noise.clusters)
        tmp.j[cc] = tmp[[1]]
        best.match[cc] = tmp[[2]]
        
        # control. check jaccard(clusters0, best.match) == tmp.j
        prots0 = unlist(strsplit(clusters0$cluster[ii], ";"))
        prots = unlist(strsplit(best.match[cc], ";"))
        jcalc = sum(prots0 %in% prots)/(length(unique(c(prots0, prots))))
        #print(paste(jcalc, tmp.j[cc]))
        if (!jcalc==tmp.j[cc]) error
      }
    }
    
    clusters0$n.removed[ii] = paste(n.changed[,ii], collapse = ";")
    clusters0$reproducibility.J[ii] = mean(tmp.j, na.rm=T)
    clusters0$all.iterations.repJ[ii] = paste(tmp.j, collapse = ";")
    clusters0$best.match[ii] = paste(best.match, collapse = " | ")
    
    # find best50, best 75
    # nodes that are in >0.50 (>0.75) of best.matches
    allprots = unlist(strsplit(unlist(strsplit(clusters0$best.match[ii], " | ", fixed=T)), ";"))
    freq = table(allprots)
    freq50 = freq>(round(iters * 0.50))
    freq75 = freq>(round(iters * 0.75))
    freq90 = freq>(round(iters * 0.90))
    if (sum(freq50)>1) clusters0$best50[ii] = paste(sort(names(freq)[freq50]), collapse = ";")
    if (sum(freq75)>1) clusters0$best75[ii] = paste(sort(names(freq)[freq75]), collapse = ";")
    if (sum(freq90)>1) clusters0$best90[ii] = paste(sort(names(freq)[freq90]), collapse = ";")
  }
  
  return(clusters0)
}