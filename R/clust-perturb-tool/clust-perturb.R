
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
    cc = 0
    for (jj in 1:iters) {
      for (kk in 1:length(noise)) {
        noise.clusters = clusters.noise$cluster[clusters.noise$iter==jj & 
                                                  clusters.noise$noise_mag==noise[kk]]
        cc = cc+1
        tmp.j[cc] = calcA(clusters0$cluster[ii], noise.clusters)
      }
    }
    
    clusters0$reproducibility.J[ii] = mean(tmp.j, na.rm=T)
  }
  
  return(clusters0)
}

