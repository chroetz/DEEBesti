derivFunGaussianProcess <- function(u, parms, bandwidth, regulation) {
  knn <- parms$knnFun(u)
  state <- parms$trajs$state[knn$idx, , drop=FALSE]
  deriv <- parms$trajs$deriv[knn$idx, , drop=FALSE]
  kernelMatrix <- expKernelMatrix(state, bandwidth, regulation)
  kernelVector <- expKernelVector(knn$distSqr, bandwidth)
  crossprod(kernelVector, solve.default(kernelMatrix, deriv))
}

