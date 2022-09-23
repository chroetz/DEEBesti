buildDerivFun <- function(name) {
  switch(
    name,
    NearestNeighbor = derivFunNearestNeighbor,
    InterpolKNN = derivFunInterpolKNN,
    KernelKNN = derivFunKernelKNN,
    LocalConst = derivFunLocalConst,
    LocalLinear = derivFunLocalLinear,
    stop("Unknown derivFun name: ", name)
  )
}

selectDerivFunHyperParms <- function(hyperParms) {
  derivFunParms <- list()
  if (hyperParms$derivFun == "NearestNeighbor") {
  } else if (hyperParms$derivFun == "InterpolKNN") {
    derivFunParms$k <- hyperParms$derivFunK
    derivFunParms$p <- hyperParms$derivFunP
  } else if (hyperParms$derivFun == "KernelKNN") {
    derivFunParms$k <- hyperParms$derivFunK
    derivFunParms$bw <- hyperParms$derivFunBw
    derivFunParms$kernel <- getKernel(hyperParms$derivFunKernel)
  } else if (hyperParms$derivFun %in% c("LocalConst", "LocalLinear")) {
    derivFunParms$bw <- hyperParms$derivFunBw
    derivFunParms$kernel <- getKernel(hyperParms$derivFunKernel)
  } else {
    stop("Unknown derivFun name ", hyperParms$derivFun)
  }
  derivFunParms
}


derivFunNearestNeighbor <- function(t, u, parms) {
  dstSqr <- rowSums((parms$trajs$state - rep(u, each=nrow(parms$trajs$state)))^2)
  du <- parms$trajs$deriv[which.min(dstSqr[-length(dstSqr)]), ]
  list(du)
}


derivFunInterpolKNN <- function(t, u, parms) {
  dstSqr <- rowSums((parms$trajs$state - rep(u, each=nrow(parms$trajs$state)))^2)
  sel <- rank(dstSqr) <= parms$derivFun$k
  dus <- parms$trajs$deriv[sel, , drop=FALSE]
  if (parms$derivFun$p > 0) { # p > 0 (p==2): interpolation weights, for p = 0: no weights
    dstSel <- dstSqr[sel]^(parms$derivFun$p/2)
    w <- 1/dstSel
    w <- w / sum(w)
    w[!is.finite(w)] <- 1/length(w)
    w <- w / sum(w)
    du <- colMeans(dus * w)
  } else {
    du <- colMeans(dus)
  }
  list(du)
}

derivFunKernelKNN <- function(t, u, parms) {
  dstSqr <- rowSums((parms$trajs$state - rep(u, each=nrow(parms$trajs$state)))^2)
  sel <- rank(dstSqr, ties.method="first") <= parms$derivFun$k
  dus <- parms$trajs$deriv[sel, , drop=FALSE]
  dstSqrSel <- dstSqr[sel]
  w <- parms$derivFun$kernel(sqrt(dstSqrSel) / parms$derivFun$bw)
  w <- w / sum(w)
  w[!is.finite(w)] <- 1/length(w)
  w <- w / sum(w)
  du <- colMeans(dus * w)
  list(du)
}

derivFunLocalConst <- function(t, u, parms) {
  dst <- sqrt(rowSums((rep(u, each=nrow(parms$trajs$state)) - parms$trajs$state)^2))
  w <- parms$derivFun$kernel(dst / parms$derivFun$bw)
  w <- w / sum(w)
  w[!is.finite(w)] <- 1/length(w)
  w <- w / sum(w)
  du <- colSums(parms$trajs$deriv * w)
  du[is.na(du)] <- 0
  list(du)
}

derivFunLocalLinear <- function(t, u, parms) {
  # TODO
  stop("derivFunLocalLinear is not implemented yet")
}

