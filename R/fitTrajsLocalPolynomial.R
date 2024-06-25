fitTrajsLocalPolynomial <- function(obs, opts) {
  opts <- asOpts(opts, c("LocalPolynomial", "FitTrajs"))
  outTime <- interpolateTimeAdaptive(obs, opts$interSteps)
  mapTrajs2Trajs(
    obs,
    fitTrajLocalPolynomial,
    outTime = outTime,
    bandwidth = opts$bandwidth,
    kernel = getKernel(opts$kernel),
    degree = opts$degree)
}


fitTrajLocalPolynomial <- function(traj, outTime, bandwidth, kernel, degree) {
  X <- outer(traj$time, 0:degree, `^`)
  m <- length(outTime)
  estiState <- matrix(nrow = m, ncol = ncol(traj$state))
  estiDeriv <- matrix(nrow = m, ncol = ncol(traj$state))
  for (j in seq_len(m)) {
    tm <- outTime[j]
    dist <- abs(traj$time-tm)
    # ignore all times that are too far away
    sel <- dist / bandwidth < 10 # TODO: this should be an option, or depend on the kernel (its support)
    w <- kernel(dist[sel] / bandwidth)
    Xsel <- X[sel, , drop=FALSE]
    Xw <- Xsel * w
    beta <- DEEButil::saveSolve(crossprod(Xw, Xsel), crossprod(Xw, traj$state[sel, , drop=FALSE]))
    psiState <- tm^(0:degree)
    psiDeriv <- cbind(0, tm^(0:(degree-1))*(1:degree))
    estiState[j, ] <- psiState %*% beta
    estiDeriv[j, ] <- psiDeriv %*% beta
  }
  makeTrajs(
    time = outTime,
    state = estiState,
    deriv = estiDeriv)
}
