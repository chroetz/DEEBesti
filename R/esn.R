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

createEsn <- function(size, inDim, degree, spectralRadius, inWeightScale, bias, seed) {

  set.seed(seed)

  inWeightMatrix <- matrix(
      inWeightScale * stats::rnorm(size * (inDim + 1)),
      nrow = size, ncol = inDim + 1)
  if (degree == 0 || spectralRadius == 0) {
    reservoirWeightMatrix <- matrix(0, size, size)
  } else {
    tmpWeightMatrix <- sampleWeightMatrix(size, deg = degree)
    reservoirWeightMatrix <-
      tmpWeightMatrix * spectralRadius / spectralRadius(tmpWeightMatrix)
  }

  return(lst(
    inWeightMatrix,
    reservoirWeightMatrix,
    size,
    inDim,
    bias))
}

trainEsn <- function(esn, obs, l2Penalty, warmUpLen, initReservoirScale, timeStepAsInput, skip = 0) {

  if (skip > 0) {
    m <- length(getTrajIds(obs))
    stopifnot(1:m == getTrajIds(obs))
    obs <-
      obs |>
      dplyr::arrange(trajId, time) |>
      dplyr::mutate(trajId = trajId + (0:skip)*m) |>
      dplyr::arrange(trajId, time)
  }

  reservoirSeries <- mapTrajs2Trajs(obs, \(traj) {

    if (timeStepAsInput) {
      timeSteps <- diff(traj$time)
      timeSteps <- c(timeSteps, mean(timeSteps))
    }

    trajReservoirSeries <- matrix(NA_real_, nrow = nrow(traj), ncol = esn$size)
    reservoir <- stats::rnorm(esn$size, sd = initReservoirScale/esn$size)

    for (i in seq_len(nrow(traj))) {
      if (timeStepAsInput) {
        v <- c(esn$bias, traj$state[i, ], timeSteps[i])
      } else {
        v <- c(esn$bias, traj$state[i, ])
      }
      reservoir <- tanh(
        esn$inWeightMatrix %*% v +
        esn$reservoirWeightMatrix %*% reservoir)
      trajReservoirSeries[i,] <- reservoir
    }

    out <- makeTrajs(
      time = traj$time,
      state = trajReservoirSeries)

    return(out)
  })

  regressionOut <- do.call(rbind, applyTrajId(obs, \(traj) traj$state[-1,, drop=FALSE]))
  regressionIn <- do.call(rbind, applyTrajId(reservoirSeries, \(traj) traj$state[-nrow(traj),, drop=FALSE]))

  X <- cbind(1, regressionIn)
  XTX <- crossprod(X)
  diag(XTX) <- diag(XTX) + c(0, rep(l2Penalty, esn$size))

  outWeightMatrix <- DEEButil::saveSolve(
    XTX,
    crossprod(X, regressionOut))

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE) # mean timeStep

  return(c(esn, lst(
    outWeightMatrix,
    timeStep,
    timeStepAsInput,
    reservoir = reservoirSeries$state,
    states = obs$state)))
}

predictEsn <- function(esn, startState, len = NULL, startTime = 0, timeRange = NULL) {

  if (is.null(timeRange)) {
    stopifnot(length(len) == 1, len >= 0)
    time <- startTime + (0:len)*esn$timeStep
  } else {
    stopifnot(length(timeRange) == 2)
    time <- seq(timeRange[1], timeRange[2], by = esn$timeStep)
    if (time[length(time)] < timeRange[2]) {
      time <- c(time, time[length(time)] + esn$timeStep)
    }
    len <- length(time) - 1
  }

  outStates <- matrix(NA_real_, nrow = len+1, ncol = ncol(esn$outWeightMatrix))
  outStates[1, ] <- startState

  # Decide how to initialize the reservoir
  iStart <- DEEButil::whichMinDist(esn$states, startState)
  if (sum((esn$states[iStart,] - startState)^2) < sqrt(.Machine$double.eps)) {
    startReservoir <- esn$reservoir[iStart, ]
    reservoir <- startReservoir
  } else {
    reservoir <- tanh(esn$inWeightMatrix %*% c(esn$bias, startState))
  }

  for (i in seq_len(len)) {
    x <- crossprod(esn$outWeightMatrix, c(1, reservoir))
    outStates[i+1,] <- x
    if (esn$timeStepAsInput) {
      v <- c(esn$bias, x, esn$timeStep)
    } else {
      v <- c(esn$bias, x)
    }
    reservoir <- tanh(
      esn$inWeightMatrix %*% v +
      esn$reservoirWeightMatrix %*% reservoir)
  }

  outTrajs <- makeTrajs(
    time = time,
    state = outStates)

  return(outTrajs)
}


predictEsnDeriv <- function(esn, states, derivOrder) {
  t(apply(states, 1, \(s) {
    predictedStates <- predictEsn(esn, s, len = derivOrder)$state
    polyInterpCoeffs <- polynomialInterpolation(esn$timeStep * 0:derivOrder, predictedStates)
    polyInterpCoeffs[2,] # derivative at 0 of polynomial is linear coefficient (second coeff)
  }))
}

polynomialInterpolation <- function(x, y) {
  stopifnot(length(x) == nrow(y))
  p <- length(x) - 1
  X <- outer(x, (0:p), `^`)
  coeff <- DEEButil::saveSolve(crossprod(X), crossprod(X, y))
  return(coeff)
}
