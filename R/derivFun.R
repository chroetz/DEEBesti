derivFunNearestNeighbor <- function(t, u, parms) {
  dst <- rowSums((parms$state - rep(u, each=nrow(parms$state)))^2)
  du <- parms$deriv[which.min(dst[-length(dst)]), ]
  list(du)
}
