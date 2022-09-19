derivFunNearestNeighbor <- function(t, u, parms) {
  dstSqr <- rowSums((parms$state - rep(u, each=nrow(parms$state)))^2)
  du <- parms$deriv[which.min(dstSqr[-length(dstSqr)]), ]
  list(du)
}


derivFunLocalConst <- function(t, u, parms) {
  dst <- sqrt(rowSums((rep(u, each=nrow(parms$state)) - parms$state)^2))
  weights <- parms$kernel(dst / parms$bw)
  weights <- weights/sum(weights)
  du <- colSums(parms$deriv * weights)
  du[is.na(du)] <- 0
  list(du)
}


buildDerivFun <- function(name) {
  switch(
    name,
    NearestNeighbor = derivFunNearestNeighbor,
    LocalConst = derivFunLocalConst,
    stop("Unknown derivFun name: ", name)
  )
}
