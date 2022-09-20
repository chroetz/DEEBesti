buildDerivFun <- function(name) {
  switch(
    name,
    NearestNeighbor = derivFunNearestNeighbor,
    WeightedKNN = derivFunWeightedKNN,
    LocalConst = derivFunLocalConst,
    LocalLinear = derivFunLocalLinear,
    stop("Unknown derivFun name: ", name)
  )
}


derivFunNearestNeighbor <- function(t, u, parms) {
  dstSqr <- rowSums((parms$trajs$state - rep(u, each=nrow(parms$trajs$state)))^2)
  du <- parms$trajs$deriv[which.min(dstSqr[-length(dstSqr)]), ]
  list(du)
}


derivFunWeightedKNN <- function(t, u, parms) {
  dstSqr <- rowSums((parms$trajs$state - rep(u, each=nrow(parms$trajs$state)))^2)
  sel <- rank(dstSqr) <= parms$derivFun$k
  dus <- parms$trajs$deriv[sel, , drop=FALSE]
  if (parms$derivFun$p > 0) { # p > 0 (p==2): interpolation weights, for p = 0: no weights
    dstSqrSel <- dstSqr[sel]^(parms$derivFun$p/2)
    w <- 1/dstSqrSel
    w <- w / sum(w)
    w[is.na(w)] <- 1
    du <- colMeans(rep(w, time=ncol(dus)) * dus)
  } else {
    du <- colMeans(dus)
  }
  list(du)
}

derivFunLocalConst <- function(t, u, parms) {
  dst <- sqrt(rowSums((rep(u, each=nrow(parms$trajs$state)) - parms$trajs$state)^2))
  weights <- parms$derivFun$kernel(dst / parms$derivFun$bw)
  weights <- weights/sum(weights)
  du <- colSums(parms$trajs$deriv * weights)
  du[is.na(du)] <- 0
  list(du)
}

derivFunLocalLinear <- function(t, u, parms) {
  # TODO
  stop("derivFunLocalLinear is not implemented yet")
}

