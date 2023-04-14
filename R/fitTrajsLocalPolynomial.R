fitTrajsLocalPolynomial <- function(obs, opts) {
  opts <- asOpts(opts, c("LocalPolynomial", "FitTrajs"))
  outTime <- interpolateTime(obs$time, opts$interSteps)
  mapTrajs2Trajs(
    obs,
    fitTrajLocalPolynomial,
    outTime = outTime,
    bandwidth = opts$bandwidth,
    kernel = getKernel(opts$kernel),
    degree = opts$degree)
}


fitTrajLocalPolynomial <- function(traj, outTime, bandwidth, kernel, degree) {
  m <- length(outTime)
  X <- outer(traj$time, 0:degree, `^`)
  Ws <- kernel(outer(traj$time, outTime, `-`) / bandwidth)
  XTWs <- array(
    rep(t(X), times = m) * rep(Ws, each = degree + 1),
    dim = c(degree + 1, length(traj$time), m))
  psiState <- outer(outTime, 0:degree, `^`)
  psiDeriv <- cbind(0, outer(outTime, 0:(degree-1), \(u, s) u^s*(s+1)))
  estiState <- matrix(nrow = m, ncol = ncol(traj$state))
  estiDeriv <- matrix(nrow = m, ncol = ncol(traj$state))
  regu <- .Machine$double.eps
  for (j in 1:m) {
    B <- XTWs[,,j] %*% X
    a <- XTWs[,,j] %*% traj$state
    Z <- DEEButil::saveSolve(B, a)
    estiState[j, ] <- psiState[j,] %*% Z
    estiDeriv[j, ] <- psiDeriv[j,] %*% Z
  }
  makeTrajs(
    time = outTime,
    state = estiState,
    deriv = estiDeriv)
}
