derivFunGaussianProcess <- function(
    u, parms, bandwidth = 1, regulation = sqrt(.Machine$double.eps),
    kMax = 100
) {
  # TODO test speed loss
  idx <- lookAround(parms$look, u, kMax)
  stopifnot(length(idx) <= kMax)
  state <- parms$trajs$state[idx, , drop=FALSE]
  deriv <- parms$trajs$deriv[idx, , drop=FALSE]
  # state <- parms$trajs$state
  # deriv <- parms$trajs$deriv
  # if (nrow(state) > kMax) {
  #   dst <- distToVec(state, u)
  #   sel <- .Internal(rank(dst, length(dst), "min")) <= kMax
  #   state <- state[sel, , drop=FALSE]
  #   deriv <- deriv[sel, , drop=FALSE]
  # }
  n <- nrow(state)
  if (n == 0) return(rep(0, length(u)))

  K <- expKernelMatrix(state, bandwidth, regulation)
  k <- expKernelVector(state, u, bandwidth)
  crossprod(k, solve.default(K, deriv))
}

