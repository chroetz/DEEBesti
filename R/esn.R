sampleWeightMatrix <- function(n, p = NULL, deg = NULL) {

  n <- as.integer(n)
  stopifnot(is.numeric(n), length(n) == 1, !is.na(n), n >= 0)
  stopifnot(is.null(p) || is.null(deg))
  if (is.null(p)) p <- deg / n
  stopifnot(is.numeric(p), length(p) == 1, p >= 0, p <= 1)

  adjacency <- stats::runif(n*n) <= p
  weights <- stats::runif(n*n, min=-1, max=1)
  weightMatrix <- matrix(adjacency*weights, nrow = n)

  return(weightMatrix)
}

spectralRadius <- function(mat) {
  max(abs(eigen(mat, symmetric=FALSE, only.values=TRUE)$values))
}

createEsn <- function(size, inDim, degree, spectralRadius, inWeightScale, bias) {
  inWeightMatrix <- matrix(
      inWeightScale * stats::rnorm(size * (inDim + 1)),
      nrow = size, ncol = inDim + 1)
  tmpWeightMatrix <- sampleWeightMatrix(size, deg = degree)
  reservoirWeightMatrix <-
    tmpWeightMatrix * spectralRadius / spectralRadius(tmpWeightMatrix)

  return(lst(
    inWeightMatrix,
    reservoirWeightMatrix,
    size,
    inDim,
    bias))
}

trainEsn <- function(esn, obs, l2Penalty, warmUpLen, initReservoirScale) {

  reservoir <- rnorm(esn$size, sd = initReservoirScale/esn$size)
  reservoirSeries <- matrix(NA_real_, nrow = nrow(obs), ncol = esn$size)
  reservoirSeries[1,] <- reservoir

  for (i in seq_len(nrow(obs)-1)) {
    reservoir <- tanh(
      esn$inWeightMatrix %*% c(esn$bias, obs$state[i, ]) +
      esn$reservoirWeightMatrix %*% reservoir)
    reservoirSeries[i+1,] <- reservoir
  }

  X <- cbind(1, reservoirSeries[-seq_len(warmUpLen),])
  XTX <- crossprod(X)
  diag(XTX) <- diag(XTX) + c(0, rep(l2Penalty, esn$size))

  outWeightMatrix <- solve.default(
    XTX,
    crossprod(X, obs$state[-seq_len(warmUpLen),, drop=FALSE]))

  return(c(esn, lst(
    outWeightMatrix,
    obs[-seq_len(warmUpLen),],
    reservoirSeries[-seq_len(warmUpLen),])))
}

predictEsn <- function(esn, initReservoir, len) {

  allX <- matrix(NA_real_, nrow = len, ncol = ncol(esn$outWeightMatrix))

  reservoir <- initReservoir
  x <- crossprod(esn$outWeightMatrix, c(1, reservoir))
  allX[1,] <- x

  for (i in seq_len(len-1)) {
    reservoir <- tanh(
      esn$inWeightMatrix %*% c(esn$bias, x) +
      esn$reservoirWeightMatrix %*% reservoir)
    x <- crossprod(esn$outWeightMatrix, c(1, reservoir))
    allX[i+1,] <- x
  }

  return(allX)
}
