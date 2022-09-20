buildDerivFun <- function(name) {
  switch(
    name,
    NearestNeighbor = derivFunNearestNeighbor,
    LocalConst = derivFunLocalConst,
    stop("Unknown derivFun name: ", name)
  )
}


derivFunNearestNeighbor <- function(t, u, parms) {
  dstSqr <- rowSums((parms$trajs$state - rep(u, each=nrow(parms$trajs$state)))^2)
  du <- parms$trajs$deriv[which.min(dstSqr[-length(dstSqr)]), ]
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

