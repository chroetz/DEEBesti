buildDerivFun <- function(opts) {
  opts <- asOpts(opts, "DerivFun")
  name <- getClassAt(opts, 2)
  derivFun <- switch(
    name,
    Null = \(t, u, parms) rep(0, length(u)),
    NearestLine = \(t, u, parms) derivFunNearestLine(
      u, parms, opts$target),
    Knn = \(t, u, parms) derivFunKnn(
      u, parms),
    GlobalLm = \(t, u, parms) derivFunGlobalLm(
      u, parms),
    LocalConst = \(t, u, parms) derivFunLocalConst(
      u, parms, bw = opts$bandwidth, kernel = getKernel(opts$kernel)),
    LocalLinear = \(t, u, parms) derivFunLocalLinear(
      u, parms, bw = opts$bandwidth, kernel = getKernel(opts$kernel)),
    GaussianProcess = \(t, u, parms) derivFunGaussianProcess(
      u, parms, bandwidth = opts$bandwidth, regulation = opts$regulation),
    InverseDistance = \(t, u, parms) derivFunInverseDistance(
      u, parms, p = opts$p),
    stop("Unknown derivFun name: ", name)
  )
  for (modOpts in opts$modifierList$list) {
    name <- getClassAt(modOpts, 2)
    derivFun <- switch(
      name,
      Attraction = addAttraction(derivFun, modOpts$attraction),
      stop("Unknown modifier name: ", name)
    )
  }
  function(t, u, parms) list(derivFun(t, u, parms))
}


derivFunNearestLine <- function(u, parms, target) {
  trajs <- parms$trajs
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
  if (any(is.na(u))) return(rep(NA, length(u)))
  knn <- parms$knnFun(u)
  if (any(knn$idx == 0)) return(rep(NA, length(u)))
  deriv <- parms$trajs$deriv[knn$idx, , drop=FALSE]
  du <- colMeans(deriv)
  return(du)
}


addAttraction <- function(derivFun, attraction) {
  force(derivFun)
  force(attraction)
  function(t, u, parms) {
    du <- derivFun(t, u, parms)
    norm <- sqrt(sum(du^2))
    knn <- parms$attactionKnnFun(u)
    nearTrajsState <- colMeans(parms$trajs$state[knn$idx, , drop=FALSE])
    du <- du + attraction * (nearTrajsState - u)
    du <- du * norm / sqrt(sum(du^2))
    return(du)
  }
}


derivFunGaussianProcess <- function(u, parms, bandwidth, regulation) {
  if (any(is.na(u))) return(rep(NA, length(u)))
  knn <- parms$knnFun(u)
  if (any(knn$idx == 0)) return(rep(NA, length(u)))
  state <- parms$trajs$state[knn$idx, , drop=FALSE]
  deriv <- parms$trajs$deriv[knn$idx, , drop=FALSE]
  kernelMatrix <- expKernelMatrix(state, bandwidth, regulation)
  kernelVector <- expKernelVectorFromDistSqr(knn$distSqr, bandwidth)
  crossprod(kernelVector, solve.default(kernelMatrix, deriv))
}


derivFunInverseDistance <- function(u, parms, p) {
  if (any(is.na(u))) return(rep(NA, length(u)))
  knn <- parms$knnFun(u)
  if (any(knn$idx == 0)) return(rep(NA, length(u)))
  deriv <- parms$trajs$deriv[knn$idx, , drop=FALSE]
  w <- 1/(knn$distSqr^(p/2))
  du <- weightedMean(deriv, w)
  return(du)
}


derivFunLocalConst <- function(u, parms, kernel, bw) {
  if (any(is.na(u))) return(rep(NA, length(u)))
  knn <- parms$knnFun(u)
  if (any(knn$idx == 0)) return(rep(NA, length(u)))
  deriv <- parms$trajs$deriv[knn$idx, , drop=FALSE]
  w <- kernel(sqrt(knn$distSqr) / bw)
  du <- weightedMean(deriv, w)
  return(du)
}


derivFunLocalLinear <- function(u, parms, kernel, bw) {
  if (any(is.na(u))) return(rep(NA, length(u)))
  knn <- parms$knnFun(u)
  if (any(knn$idx == 0)) return(rep(NA, length(u)))
  state <- parms$trajs$state[knn$idx, , drop=FALSE]
  deriv <- parms$trajs$deriv[knn$idx, , drop=FALSE]
  w <- kernel(sqrt(knn$distSqr) / bw)
  X <- cbind(1, state)
  Xw <- X * w

  regu <- .Machine$double.eps
  while(TRUE) {
    tryCatch(
      {
        beta <- solve.default(crossprod(Xw, X) + regu * diag(ncol(X)), crossprod(Xw, deriv))
        break
      },
      error = function(cond) regu <<- regu * 10
    )
    cat("*")
    if (regu > sqrt(sqrt(.Machine$double.eps))) return(rep(NA, length(u)))
  }

  crossprod(c(1, u), beta)
}


weightedMean <- function(x, w) {
  w <- w / sum(w)
  w[!is.finite(w)] <- 1/length(w)
  w <- w / sum(w)
  colSums(x * w)
}

