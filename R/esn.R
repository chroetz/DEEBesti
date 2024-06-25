sampleWeightMatrix <- function(n, p = NULL, deg = NULL) {

  n <- as.integer(n)
  stopifnot(is.numeric(n), length(n) == 1, !is.na(n), n >= 0)
  stopifnot(is.null(p) || is.null(deg))
  if (is.null(p)) p <- deg / n
  stopifnot(is.numeric(p), length(p) == 1, p >= 0, p <= 1)

  adjacency <- stats::runif(n*n) <= p
  weights <- stats::runif(n*n, min=-1, max=1)
  weightMatrix <- matrix(adjacency*weights, nrow = n)

  return(weightMatrix)
}


spectralRadius <- function(mat) {
  max(abs(eigen(mat, symmetric=FALSE, only.values=TRUE)$values))
}


createEsn <- function(opts, inDim) {

  opts <- asOpts(opts, c("Esn", "Propagator", "HyperParms"))

  set.seed(opts$seed)

  inWeightMatrix <- matrix(
      opts$inWeightScale * stats::rnorm(opts$size * (inDim + 1)),
      nrow = opts$size, ncol = inDim + 1)
  if (opts$degree == 0 || opts$spectralRadius == 0) {
    reservoirWeightMatrix <- matrix(0, opts$size, opts$size)
  } else {
    tmpWeightMatrix <- sampleWeightMatrix(opts$size, deg = opts$degree)
    reservoirWeightMatrix <-
      tmpWeightMatrix * opts$spectralRadius / spectralRadius(tmpWeightMatrix)
  }

  return(lst(
    inWeightMatrix,
    reservoirWeightMatrix,
    inDim))
}


trainEsn <- function(parms, obs, opts) {

  opts <- asOpts(opts, c("Esn", "Propagator", "HyperParms"))

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

    trajReservoirSeries <- matrix(NA_real_, nrow = nrow(traj), ncol = opts$size)
    reservoir <- stats::rnorm(opts$size, sd = opts$initReservoirScale/opts$size)

    for (i in seq_len(nrow(traj))) {
      if (opts$timeStepAsInput) {
        v <- c(opts$bias, traj$state[i, ], timeSteps[i])
      } else {
        v <- c(opts$bias, traj$state[i, ])
      }
      reservoir <- tanh(
        parms$inWeightMatrix %*% v +
        parms$reservoirWeightMatrix %*% reservoir)
      trajReservoirSeries[i,] <- reservoir
    }

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

  outWeightMatrix <- DEEButil::saveSolve(
    XTX,
    crossprod(X, regressionOut))

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE) # mean timeStep

  return(c(parms, lst(
    outWeightMatrix,
    timeStep,
    reservoir = reservoirSeries$state,
    states = obs$state)))
}


predictEsn <- function(parms, opts, startState, len) {

  opts <- asOpts(opts, c("Esn", "Propagator", "HyperParms"))

  outStates <- matrix(NA_real_, nrow = len+1, ncol = length(startState))
  outStates[1, ] <- startState

  if (!hasValue(parms)) return(outStates)

  # Decide how to initialize the reservoir
  iStart <- DEEButil::whichMinDist(parms$states, startState)
  if (sum((parms$states[iStart,] - startState)^2) < sqrt(.Machine$double.eps)) {
    cat("Found startState in training data. Use it to initialize Reservoir.\n")
    startReservoir <- parms$reservoir[iStart, ]
    reservoir <- startReservoir
  } else {
    cat("Did not find startState in training data. Use 0 reservoir.\n")
    if (opts$timeStepAsInput) {
      v <- c(opts$bias, startState, parms$timeStep)
    } else {
      v <- c(opts$bias, startState)
    }
    reservoir <- tanh(parms$inWeightMatrix %*% v)
  }

  prevState <- startState
  for (i in seq_len(len)) {
    prediction <- as.vector(c(1, reservoir) %*% parms$outWeightMatrix)
    newState <- getPropagatorNextState(prevState, parms$timeStep, prediction, opts$targetType)
    outStates[i+1,] <- newState
    if (opts$timeStepAsInput) {
      v <- c(opts$bias, newState, parms$timeStep)
    } else {
      v <- c(opts$bias, newState)
    }
    reservoir <- tanh(
      parms$inWeightMatrix %*% v +
      parms$reservoirWeightMatrix %*% reservoir)
    prevState <- newState
  }

  return(outStates)
}
