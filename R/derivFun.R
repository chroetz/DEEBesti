buildDerivFun <- function(opts) {
  opts <- asOpts(opts, "DerivFun")
  name <- getClassAt(opts, 2)
  derivFunUnlisted <- switch(
    name,
    Null = \(t, u, parms) rep(0, length(u)),
    NearestNeighbor = \(t, u, parms) derivFunNearestNeighbor(
      u, parms),
    InterpolKNN = \(t, u, parms) derivFunInterpolKNN(
      u, parms, p = opts$p, k = opts$k),
    KernelKNN = \(t, u, parms) derivFunKernelKNN(
      u, parms, k = opts$k, bw = opts$bandwidth, kernel = getKernel(opts$kernel)),
    LocalConst = \(t, u, parms) derivFunLocalConst(
      u, parms, bw = opts$bandwidth, kernel = getKernel(opts$kernel)),
    LocalLinear = \(t, u, parms) derivFunLocalLinear(
      u, parms, bw = opts$bandwidth, kernel = getKernel(opts$kernel)),
    stop("Unknown derivFun name: ", name)
  )
  function(t, u, parms) list(derivFunUnlisted(t, u, parms))
}


derivFunNearestNeighbor <- function(u, trajs) {
  trajs$deriv[whichMinDist(trajs$state, u), ]
}


derivFunInterpolKNN <- function(u, trajs, p, k) { # BEWARE: this will introduce discontinuities
  dst <- distToVec(trajs$state, u)
  sel <- rank(dst) <= k
  dus <- trajs$deriv[sel, , drop=FALSE]
  if (p > 0) { # p > 0 (p==2): interpolation weights, for p = 0: no weights
    dstSel <- dst[sel]^p
    w <- 1/dstSel
    w <- w / sum(w)
    w[!is.finite(w)] <- 1/length(w)
    w <- w / sum(w)
    du <- colMeans(dus * w)
  } else {
    du <- colMeans(dus)
  }
  return(du)
}

derivFunKernelKNN <- function(u, trajs, k, kernel, bw) {
  dst <- distToVec(trajs$state, u)
  sel <- rank(dst, ties.method="first") <= k
  dus <- trajs$deriv[sel, , drop=FALSE]
  dstSel <- dst[sel]
  w <- kernel(dstSel / bw)
  w <- w / sum(w)
  w[!is.finite(w)] <- 1/length(w)
  w <- w / sum(w)
  du <- colMeans(dus * w)
  return(du)
}

derivFunLocalConst <- function(u, trajs, kernel, bw) {
  dst <- distToVec(trajs$state, u)
  w <- kernel(dst / bw)
  w <- w / sum(w)
  w[!is.finite(w)] <- 1/length(w)
  w <- w / sum(w)
  du <- colSums(trajs$deriv * w)
  du[is.na(du)] <- 0
  return(du)
}

derivFunLocalLinear <- function(u, trajs, kernel, bw) {
  # TODO
  stop("derivFunLocalLinear is not implemented yet")
}

