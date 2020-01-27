mymcl = function (m, infl, iter = 1000, remove.self.loops = FALSE, prune = FALSE, 
                  thresh = 1e-06, pruning.prob = 1e-06, use.sparse = NULL, 
                  verbose = FALSE) 
{
  if (nrow(m) != ncol(m)) {
    stop("input matrix must be a square matrix")
  }
  if (remove.self.loops) {
    diag(m) = 0
  }
  n = nrow(m)^2
  m <- sweep(m, 2, colSums(m), `/`)
  m[is.na(m)] <- 0
  if (prune) {
    m[m < pruning.prob] = 0
  }
  force.sparse = FALSE
  if (length(use.sparse) == 0) {
    force.sparse = ((sum(m == 0)/n) >= 0.5)
    use.sparse = force.sparse
  }
  if (use.sparse || force.sparse) {
    {
      m = Matrix(m)
      if (verbose) {
        print("sparse matrices will be used")
      }
    }
  }
  m0 <- m
  m <- m %*% m
  m <- m^infl
  m <- sweep(m, 2, colSums(m), `/`)
  m[is.na(m)] <- 0
  i = 1
  if (sum(m0 - m) != 0) {
    for (i in 2:iter) {
      m <- m %*% m
      m <- m^infl
      m <- sweep(m, 2, colSums(m), `/`)
      m[is.na(m)] <- 0
      if ((sum(m > 0 & m < 1) == 0) || (sqrt((sum((m - 
                                                   m0)^2)/n)) < thresh)) {
        break
      }
      if (prune) {
        m[m < pruning.prob] <- 0
      }
      m0 = m
    }
  }
  if (verbose) {
    print(paste("mcl converged after", i, "iterations"))
  }
  if (class(matrix) != "matrix") {
    m = as.matrix(m)
  }
  nrow <- nrow(m)
  ncol <- ncol(m)
  clusters <- vector(length = nrow, mode = "numeric")
  csums = colSums(m)
  lonely = which(csums == 0)
  clustered = which(csums > 0)
  clusters[lonely] = lonely
  attractors = sort(which(rowSums(m) > 0))
  j = 1
  lcc = length(clustered)
  unlabelled = lcc
  while (unlabelled > 0) {
    i = attractors[j]
    if (clusters[i] == 0) {
      attracts <- which(m[i, 1:ncol] > 0)
      clusters[attracts] <- i
    }
    unlabelled = sum(clusters[clustered] == 0)
    j = j + 1
    # fix bug where j > length(attractors)
    if (j > length(attractors)) unlabelled = 0
  }
  return(clusters)
}
