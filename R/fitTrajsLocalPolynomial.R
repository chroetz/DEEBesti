fitTrajsLocalPolynomial <- function(obs, opts) {
  opts <- asOpts(opts, c("LocalPolynomial", "FitTrajs"))
  outTime <- interpolateTimeAdaptive(obs, opts$interSteps)
  mapTrajs2Trajs(
    obs,
    fitTrajLocalPolynomial,
    outTime = outTime,
    bandwidth = opts$bandwidth,
    kernel = getKernel(opts$kernel),
    degree = opts$degree,
    maxNeighbors = opts$maxNeighbors)
}


fitTrajLocalPolynomial <- function(traj, outTime, bandwidth, kernel, degree, maxNeighbors) {

  stopifnot(all(diff(rank(outTime)) == 1))
  stopifnot(outTime[1] <= traj$time[1] + sqrt(.Machine$double.eps))
  timeStepTraj <- DEEBtrajs::getTimeStep(traj$time)
  timeStepOut <- DEEBtrajs::getTimeStep(outTime)
  stopifnot(timeStepOut <= timeStepTraj)

  X <- outer(traj$time, 0:degree, `^`)
  m <- length(outTime)
  n <- nrow(traj)
  estiState <- matrix(nrow = m, ncol = ncol(traj$state))
  estiDeriv <- matrix(nrow = m, ncol = ncol(traj$state))
  i <- 1
  s <- floor(maxNeighbors/2)
  sel <- c(rep(TRUE, s), rep(FALSE, n - s))

  for (j in seq_len(m)) {
    tOut <- outTime[j]
    w <- kernel(abs(traj$time[sel]-tOut) / bandwidth)
    Xsel <- X[sel, , drop=FALSE]
    Xw <- Xsel * w
    beta <- DEEButil::saveSolve(crossprod(Xw, Xsel), crossprod(Xw, traj$state[sel, , drop=FALSE]))
    psiState <- tOut^(0:degree)
    psiDeriv <- cbind(0, tOut^(0:(degree-1))*(1:degree))
    estiState[j, ] <- psiState %*% beta
    estiDeriv[j, ] <- psiDeriv %*% beta

    if (tOut > traj$time[i]) {
      i <- i+1
      if (i+s-1 <= n) sel[i+s-1] <- TRUE
      if (i-s-1 >= 1) sel[i-s-1] <- FALSE
    }
  }

  makeTrajs(
    time = outTime,
    state = estiState,
    deriv = estiDeriv)
}


