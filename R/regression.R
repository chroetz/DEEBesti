trainPropagatorRegression <- function(x, y, opts)  {
  knnFun <- FastKNN::buildKnnFunction(x, opts$neighbors)
  return(lst(x, y, knnFun))
}


inferPropagatorRegression <- function(parms, opts, xout) {
  fit <- matrix(NA_real_, nrow=nrow(xout), ncol=ncol(parms$y))
  for (i in seq_len(nrow(xout))) {
    xouti <- xout[i,]
    knn <- parms$knnFun(xouti)
    x <- parms$x[knn$idx, , drop=FALSE]
    y <- parms$y[knn$idx, , drop=FALSE]
    fit[i, ] <- predictPropagatorRegression(knn$distSqr, x, y, opts)
  }
  return(fit)
}


predictPropagatorRegression <- function(distSqr, x, y, xout, opts) {
  name <- ConfigOpts::getClassAt(opts, 2)
  switch(
    name,
    LocalConst = predictPropagatorRegressionLocalConst(distSqr, y, opts),
    LocalLinear = predictPropagatorRegressionLocalLinear(distSqr, x, y, xout, opts),
    GaussianProcess = predictPropagatorRegressionGaussianProcess(distSqr, x, y, opts),
    stop("Unknown name", name)
  )
}


predictPropagatorRegressionLocalConst <- function(distSqr, x, y, opts) {
  w <- opts$kernel(sqrt(distSqr) / opts$bandwidth)
  yout <- weightedMean(y, w)
  return(yout)
}


predictPropagatorRegressionGaussianProcess <- function(distSqr, x, y, opts)  {
  kernelMatrix <- DEEButil::expKernelMatrix(x, opts$bandwidth, opts$regulation)
  kernelVector <- DEEButil::expKernelVectorFromDistSqr(distSqr, opts$bandwidth)
  crossprod(kernelVector, DEEButil::saveSolve(kernelMatrix, y))
}


predictPropagatorRegressionLocalLinear <- function(distSqr, x, y, xout, opts) {
  w <- opts$kernel(sqrt(distSqr) / opts$bandwidth)
  X <- cbind(1, x)
  Xw <- X * w
  beta <- DEEButil::saveSolve(crossprod(Xw, X), crossprod(Xw, y))
  crossprod(c(1, xout), beta)
}

