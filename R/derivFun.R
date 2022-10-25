buildDerivFun <- function(opts) {
  opts <- asOpts(opts, "DerivFun")
  name <- getClassAt(opts, 2)
  # KNN selection is only implemented for some derivFuns. For all others the field neighbors should be 0.
  stopifnot(opts$neighbors == 0 || name %in% c("Knn", "GaussianProcess"))
  derivFunUnlisted <- switch(
    name,
    Null = \(t, u, parms) rep(0, length(u)),
    NearestLine = \(t, u, parms) derivFunNearestLine(
      u, parms$trajs, opts$target),
    Knn = \(t, u, parms) derivFunKnn(
      u, parms),
    LocalConst = \(t, u, parms) derivFunLocalConst(
      u, parms$trajs, bw = opts$bandwidth, kernel = getKernel(opts$kernel)),
    LocalLinear = \(t, u, parms) derivFunLocalLinear(
      u, parms$trajs, bw = opts$bandwidth, kernel = getKernel(opts$kernel)),
    GaussianProcess = \(t, u, parms) derivFunGaussianProcess(
      u, parms, bandwidth = opts$bandwidth, regulation = opts$regulation),
    InverseDistance = \(t, u, parms) derivFunInverseDistance(
      u, parms$trajs, p = opts$p),
    stop("Unknown derivFun name: ", name)
  )
  function(t, u, parms) list(derivFunUnlisted(t, u, parms))
}


derivFunNearestLine <- function(u, trajs, target) {
  idxMin <- whichMinDistToPwLin(trajs$state, trajs$trajId, u)
  if (target == "pointLine") {
    # assumes deriv[idxMin,] belongs to segement idxMin;
    # interpolate for corner points
    if (idxMin %% 1 == 0 && idxMin > 1) { # corner point
      du <- (trajs$deriv[idxMin, ] + trajs$deriv[idxMin-1, ]) / 2
    } else {
      du <- trajs$deriv[floor(idxMin), ]
    }
  } else if (target == "interp") {
    # assumes deriv[idxMin,] belongs to point idxMin;
    # linearly interpolate between derivs for line
    if (idxMin %% 1 == 0) { # corner point
      du <- trajs$deriv[idxMin, ]
    } else {
      t <- idxMin - floor(idxMin)
      du <- (1-t)*trajs$deriv[floor(idxMin), ] + t*trajs$deriv[ceiling(idxMin), ]
    }
  } else if (target == "line") {
    du <- trajs$deriv[floor(idxMin), ]
  } else {
    stop("Unknown target ", target)
  }
  return(du)
}

derivFunKnn <- function(u, parms) {
  knn <- parms$knnFun(u)
  deriv <- parms$trajs$deriv[knn$idx, , drop=FALSE]
  du <- colMeans(deriv)
  return(du)
}

derivFunGaussianProcess <- function(u, parms, bandwidth, regulation) {
  knn <- parms$knnFun(u)
  state <- parms$trajs$state[knn$idx, , drop=FALSE]
  deriv <- parms$trajs$deriv[knn$idx, , drop=FALSE]
  kernelMatrix <- expKernelMatrix(state, bandwidth, regulation)
  kernelVector <- expKernelVector(knn$distSqr, bandwidth)
  crossprod(kernelVector, solve.default(kernelMatrix, deriv))
}

weightedMean <- function(x, w) {
  w <- w / sum(w)
  w[!is.finite(w)] <- 1/length(w)
  w <- w / sum(w)
  colSums(x * w)
}

derivFunInverseDistance <- function(u, trajs, p) {
  dst <- distToVec(trajs$state, u)
  w <- 1/(exp(dst*p)-1)
  du <- weightedMean(trajs$deriv, w)
  return(du)
}

derivFunLocalConst <- function(u, trajs, kernel, bw) {
  dst <- distToVec(trajs$state, u)
  w <- kernel(dst / bw)
  du <- weightedMean(trajs$deriv, w)
  return(du)
}

derivFunLocalLinear <- function(u, trajs, kernel, bw) {
  # TODO
  stop("derivFunLocalLinear is not implemented yet")
}

