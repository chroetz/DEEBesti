#' Fit local linear estimator.
#'
#' Take predictor-response pairs and estimate new responses for the same
#' predictors via a local linear estimator. Both, predictors and responses, may
#' be multidimensional.
#'
#' @param x \code{numeric(n)} or \code{matrix(, n, p)}. n observations of
#'   p-dimensional predictors.
#' @param y \code{numeric(n)} or \code{matrix(, n, k)}. n observations of
#'   k-dimensional responses.
#' @param bw \code{numeric(1)}. The bandwidth to use.
#' @param kernel \code{function}. The kernel. Maps the Euclidean distance of
#'   points in the predictor space to a weight.
#' @return A numeric n by k matrix. The estimates of response variable at the
#'   locations of the observed predictors.
#' @export
fitterLocalLinear <- function(x, y, bw, kernel) {
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  stopifnot(NROW(x) == NROW(y))
  stopifnot(NROW(x) >= 1)
  stopifnot(is.numeric(bw) && length(bw) == 1)
  stopifnot(is.function(kernel))
  x <- as.matrix(x)
  y <- as.matrix(y)
  dst <- as.matrix(stats::dist(x))
  weights <- kernel(dst/bw)
  X <- cbind(1, x)
  fit <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(y))
  for (j in 1:nrow(x)) {
    B <- crossprod(X, weights[,j]*X)
    a <- crossprod(X, weights[,j]*y)
    fit[j, ] <- X[j, ] %*% solve(B, a)
  }
  return(fit)
}


fitterLocalConst <- function(x, y, bw, kernel) {
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  stopifnot(NROW(x) == NROW(y))
  stopifnot(NROW(x) >= 1)
  stopifnot(is.numeric(bw) && length(bw) == 1)
  stopifnot(is.function(kernel))
  x <- as.matrix(x)
  y <- as.matrix(y)
  dst <- as.matrix(stats::dist(x))
  weights <- kernel(dst/bw) # TODO: This seems to be quite slow for exp-kernel.
  fit <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(y))
  for (j in 1:nrow(x)) {
    fit[j, ] <- colSums(weights[,j] * y) / sum(weights[,j])
  }
  return(fit)
}


fitterGaussianProcess <- function(x, y, bandwidth, kernel, regulation, neighbors) {
  knnFun <- FastKNN::buildKnnFunction(x, neighbors)
  fit <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(y))
  for (i in seq_len(nrow(x))) {
    knn <- knnFun(x[i,])
    kernelMatrix <- DEEButil::expKernelMatrix(x[knn$idx, , drop=FALSE], bandwidth, regulation)
    kernelVector <- DEEButil::expKernelVectorFromDistSqr(knn$distSqr, bandwidth)
    fit[i, ] <- crossprod(kernelVector, solve.default(kernelMatrix, y[knn$idx, , drop=FALSE]))
  }
  return(fit)
}
