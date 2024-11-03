createRafda <- function(obs, opts) {

  opts <- asOpts(opts, c("Rafda", "Propagator", "HyperParms"))

  inDim <- getDim(obs)
  if (opts$timeStepAsInput) {
    inDim <- inDim + 1
  }
  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE)

  inWeightMatrix <- matrix(
    opts$inWeightScale * stats::rnorm(opts$size * (inDim + 1)),
      nrow = opts$size, ncol = inDim + 1)

    if (opts$skip > 0) {
      m <- length(getTrajIds(obs))
      n <- nrow(obs)
      stopifnot(1:m == getTrajIds(obs))
      obs <-
        obs |>
        dplyr::arrange(.data$trajId, .data$time) |>
        dplyr::mutate(trajId = .data$trajId + rep_len((0:opts$skip)*m, n)) |>
        dplyr::arrange(.data$trajId, .data$time)
    }

  reservoirSeries <- mapTrajs2Trajs(obs, \(traj) {

    if (opts$timeStepAsInput) {
      timeSteps <- diff(traj$time)
      timeSteps <- c(timeSteps, mean(timeSteps))
    }

    if (opts$timeStepAsInput) {
      V <- cbind(opts$bias, traj$state, timeSteps)
    } else {
      V <- cbind(opts$bias, traj$state)
    }
    trajReservoirSeries <- tanh(tcrossprod(V, inWeightMatrix))

    out <- makeTrajs(
      time = traj$time,
      state = trajReservoirSeries)

    return(out)
  })

  regressionIn <- do.call(
    rbind,
    applyTrajId(reservoirSeries, \(traj) traj$state[-nrow(traj),, drop=FALSE]))
  regressionOut <- getPropagatorRegressionOut(obs, opts$targetType)

  X <- cbind(1, regressionIn)
  XTX <- crossprod(X)
  diag(XTX) <- diag(XTX) + c(0, rep(opts$l2Penalty, opts$size))
  XTy <- crossprod(X, regressionOut)
  outWeightMatrixLr <- DEEButil::saveSolve(XTX, XTy)
  yPred <- X %*% outWeightMatrixLr
  err <- regressionOut - yPred
  covU <- stats::cov(err)
  svd <- svd(covU)
  covUhalf <- svd$u %*% diag(sqrt(svd$d)) %*% t(svd$v)
  sdW <- opts$weightSdScale * stats::sd(outWeightMatrixLr) / sqrt(nrow(outWeightMatrixLr))
  ensembleIdxs <- seq_len(opts$ensembleSize)
  d <- ncol(obs$state)
  dr <- opts$size + 1
  dz <- d + d*dr

  # TODO apply trajs
  stopifnot(length(unique(obs$trajId)) == 1)
  if (opts$timeStepAsInput) {
    timeSteps <- diff(traj$time)
    timeSteps <- c(timeSteps, mean(timeSteps))
  }
  traj <- obs
  ua <- traj$state[1, ]
  wa <- outWeightMatrixLr
  uaEns <- lapply(ensembleIdxs, \(j) ua + covUhalf %*% stats::rnorm(d))
  waEns <- lapply(ensembleIdxs, \(j) wa + stats::rnorm(length(wa), sd = sdW))
  zf <- matrix(NA_real_, ncol=opts$ensembleSize, nrow=dz)
  for (k in seq_len(nrow(traj)-1)) {
    for (j in ensembleIdxs) {
      if (opts$timeStepAsInput) {
        v <- c(opts$bias, uaEns[[j]], timeSteps[k])
      } else {
        v <- c(opts$bias, uaEns[[j]])
      }
      features <- tanh(inWeightMatrix %*% v)
      uf <- crossprod(waEns[[j]], c(1, features))
      zf[,j] <- c(as.vector(uf), as.vector(waEns[[j]]))
    }
    zfHat <- zf - rep(rowMeans(zf), times=opts$ensembleSize)
    pt3 <- proc.time()
    Pf <- 1/(opts$ensembleSize-1)*tcrossprod(zfHat)
    pt4 <- proc.time()
    u <- traj$state[k+1, ]
    uPerturbed <-
      rep(u, times=opts$ensembleSize) +
      covUhalf %*% matrix(stats::rnorm(d*opts$ensembleSize), nrow=d)
    za <- zf + Pf[, 1:d] %*% DEEButil::saveSolve(Pf[1:d, 1:d] + covU, uPerturbed - zf[1:d, ])
    waEns <- lapply(ensembleIdxs, \(j) matrix(za[(d+1):nrow(za), j], nrow=dr, ncol=d))
    uaEns <- lapply(ensembleIdxs, \(j) za[1:d, j])
  }
  wRafda <- rowMeans(za[(d+1):nrow(za), ])
  outWeightMatrixRafda <- matrix(wRafda, nrow=dr, ncol=d)

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE) # mean timeStep

  return(c(lst(
    inWeightMatrix,
    outWeightMatrix = outWeightMatrixRafda,
    timeStep,
    states = obs$state)))
}



predictRafda <- function(parms, opts, startState, len) {

  opts <- asOpts(opts, c("Rafda", "Propagator", "HyperParms"))

  stateDim <- length(startState)

  outStates <- matrix(NA_real_, nrow = len+1, ncol = stateDim)
  outStates[1, ] <- startState
  if (any(is.na(startState))) {
    warning("startState has NA. Returning NA.", immediate.=TRUE)
    return(outStates)
  }

  if (opts$timeStepAsInput) {
    v <- c(opts$bias, startState, parms$timeStep)
  } else {
    v <- c(opts$bias, startState)
  }
  features <- tanh(parms$inWeightMatrix %*% v)

  prevState <- startState
  for (i in seq_len(len)) {
    prediction <- as.vector(c(1, features) %*% parms$outWeightMatrix)
    newState <- getPropagatorNextState(prevState, parms$timeStep, prediction, opts$targetType)
    outStates[i+1,] <- newState
    if (any(!is.finite(newState))) break
    if (opts$timeStepAsInput) {
      v <- c(opts$bias, newState, parms$timeStep)
    } else {
      v <- c(opts$bias, newState)
    }
    features <- tanh(parms$inWeightMatrix %*% v)
    prevState <- newState
  }

  return(outStates)
}
