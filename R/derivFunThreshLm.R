derivFunThreshLm <- function(u, parms) {
  d <- length(u)
  du <- vapply(
    seq_len(d),
    \(j) {
      degVecs <- parms$fits[[j]]$degVec
      coefs <- parms$fits[[j]]$coef
      features <- DEEButil::evaluateMonomials(matrix(u, nrow=1), degVecs)
      sum(features * coefs)
    },
    numeric(1))
  return(du)
}


prepareParmsThreshLm <- function(parms, opts) {

  opts <- asOpts(opts, c("ThreshLm", "DerivFun"))

  stopifnot(grepl("^poly\\d+$", opts$features))
  degree <- as.integer(substring(opts$features, 5))
  d <- getDim(parms$trajs)
  degVecs <- DEEButil::getMonomialExponents(d, degree)
  x <- DEEButil::evaluateMonomials(parms$trajs$state, degVecs)
  y <- parms$trajs$deriv

  beta <- linSolve(x, y)
  smallIndsPrev <- beta == 0
  for (k in 1:opts$iterations) {
    smallInds <- abs(beta) < opts$threshold
    if (all(smallInds == smallIndsPrev)) break
    beta[smallInds] <- 0
    for (j in 1:d) {
      bigInds <- !smallInds[,j]
      beta[bigInds, j] <- linSolve(x[,bigInds], y[,j])
    }
    smallIndsPrev <- smallInds
  }

  fits <- apply(
    beta,
    2,
    \(b) {
      sel <- abs(b) >= opts$threshold
      list(degVec = degVecs[sel, ], coef = b[sel])
    },
    simplify=FALSE)

  parms$fits <- fits
  return(parms)
}

linSolve <- function(x, y) {
  DEEButil::saveSolve(crossprod(x), crossprod(x, y))
}
