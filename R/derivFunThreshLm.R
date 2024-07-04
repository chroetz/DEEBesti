derivFunThreshLm <- function(u, parms) {
  d <- length(u)
  if (hasValue(parms$valid) && !parms$valid) {
    return(rep(NA_real_, d))
  }
  du <- vapply(
    seq_len(d),
    \(j) {
      if (length(parms$fits[[j]]$coef) == 0) return(0)
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

  degree <- opts$polyDeg
  d <- getDim(parms$trajs)
  nFeatures <- DEEButil::numberOfTermsInPoly(degree, d)
  if (nFeatures > 2000) {
    warning("Infeasible number of dimensions in ThreshLm. Giving trivial result.")
    return(c(parms, list(valid = FALSE)))
  }
  degVecs <- DEEButil::getMonomialExponents(d, degree)
  x <- DEEButil::evaluateMonomials(parms$trajs$state, degVecs)
  y <- parms$trajs$deriv

  beta <- linSolve(x, y, opts$l2Penalty)
  smallIndsPrev <- beta == 0
  if (opts$iterations > 0) {
    for (k in 1:opts$iterations) {
      smallInds <- abs(beta) < opts$threshold
      if (all(smallInds == smallIndsPrev)) break
      beta[smallInds] <- 0
      for (j in 1:d) {
        bigInds <- !smallInds[,j]
        if (any(bigInds)) {
          beta[bigInds, j] <- linSolve(x[, bigInds, drop=FALSE], y[,j], opts$l2Penalty)
        }
      }
      smallIndsPrev <- smallInds
    }
  }

  fits <- apply(
    beta,
    2,
    \(b) {
      sel <- abs(b) >= opts$threshold
      list(degVec = degVecs[sel, , drop=FALSE], coef = b[sel])
    },
    simplify=FALSE)

  parms$fits <- fits
  parms$valid <- TRUE
  return(parms)
}


linSolve <- function(x, y, l2Penalty = 0) {
  xtx <- crossprod(x)
  if (l2Penalty > 0) {
    diag(xtx) <- diag(xtx) + l2Penalty
  }
  DEEButil::saveSolve(xtx, crossprod(x, y))
}
