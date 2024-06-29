rowIdxProd <- function(mat, idxs) {
  n <- nrow(mat)
  if (n == 0) return(double(0))
  if (length(idxs) == 0) return(rep(1, n))
  if (length(idxs) == 1) return(mat[,idxs])
  y <- mat[,idxs[[1]]]
  for (i in idxs[-1]) {
    y <- y * mat[, i]
  }
  return(y)
}


createLinFeaturesOneTraj <- function(traj, pastSteps, skip, timeStepAsInput) {
  if (timeStepAsInput) {
    timeSteps <- diff(traj$time)
    timeSteps <- c(timeSteps, mean(timeSteps))
    state <- cbind(timeSteps, traj$state)
  } else {
    state <- traj$state
  }
  if (pastSteps == 0) return(state)
  n <- nrow(state)
  d <- ncol(state)
  featuresLin <- matrix(NA_real_, nrow = n, ncol = d*(1+pastSteps))
  featuresLin[,1:d] <- state
  for (j in seq_len(pastSteps)) {
    rowIdxs <- seq_len(nrow(state))
    rowIdxs <- rowIdxs - j*(skip+1)
    nEmpty <- sum(rowIdxs <= 0)
    rowIdxs <- rowIdxs[rowIdxs>=1]
    if (nEmpty >= n) next
    featuresLin[(nEmpty+1):n, d*j+(1:d)] <- state[rowIdxs, ]
  }
  return(featuresLin)
}


createLinFeaturesOneTrajLastTime <- function(state, pastSteps, skip) {
  n <- nrow(state)
  d <- ncol(state)
  if (pastSteps == 0) return(state[n, ])
  rowIdxs <- n - (0:pastSteps)*(skip+1)
  nEmpty <- sum(rowIdxs <= 0)
  rowIdxs <- rowIdxs[rowIdxs>=1]
  c(
    state[rowIdxs, ] |>
      t() |>
      as.vector(),
    rep(NA_real_, nEmpty*d))
}


createBaseFeatures <- function(trajs, timeStepAsInput,  pastSteps, skip) {

  featuresLinTrajs <- mapTrajs2Trajs(trajs, \(traj) {

    featuresLin <- createLinFeaturesOneTraj(traj, pastSteps, skip, timeStepAsInput)

    out <- makeTrajs(
      time = traj$time,
      state = featuresLin)

    return(out)

  })

  return(featuresLinTrajs)
}


createPolyFeaturesOneTraj <- function(featuresLin, polyDeg) {
  d <- ncol(featuresLin)
  degVecs <- DEEButil::getMonomialExponents(d, polyDeg)
  DEEButil::evaluateMonomials(featuresLin, degVecs)
}


createPolyFeaturesOneTrajOne <- function(featuresLin, polyDeg) {
  d <- length(featuresLin)
  degVecs <- DEEButil::getMonomialExponents(d, polyDeg)
  DEEButil::evaluateMonomials(matrix(featuresLin, nrow=1), degVecs) |> as.vector()
}


createPolyFeatures <- function(baseFeatures, polyDeg) {

  mapTrajs2Trajs(baseFeatures, \(traj) {

    features <- createPolyFeaturesOneTraj(traj$state, polyDeg)

    out <- makeTrajs(
      time = traj$time,
      state = features)

    return(out)

  })
}


createLinear <- function(obs, opts) {

  opts <- asOpts(opts, c("Linear", "Propagator", "HyperParms"))

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE)

  baseFeatures <- createBaseFeatures(obs, opts$timeStepAsInput,  opts$pastSteps, opts$skip)
  nFeatures <- DEEButil::numberOfTermsInPoly(opts$polyDeg, ncol(baseFeatures$state))
  if (nFeatures > 2000 || nFeatures > nrow(obs$state)) {
    warning("Infeasible number of dimensions in LinearPropagator. Giving trivial result.")
    return(lst(timeStep, obs, valid = FALSE))
  }
  featureSeries <- createPolyFeatures(baseFeatures, opts$polyDeg)

  regressionOut <- getPropagatorRegressionOut(obs, opts$targetType)
  regressionIn <- do.call(
    rbind,
    applyTrajId(featureSeries, \(traj) traj$state[-nrow(traj),, drop=FALSE]))

  naRows <- rowSums(is.na(regressionIn) > 0)
  regressionIn <- regressionIn[!naRows, ]
  regressionOut <- regressionOut[!naRows, ]

  X <- regressionIn
  XTX <- crossprod(X)
  if (opts$targetType == "deriv") {
    diag(XTX) <- diag(XTX) + rep(opts$l2Penalty, ncol(regressionIn))
  } else if (opts$targetType == "state") {
    diag(XTX) <- diag(XTX) + c(0, rep(opts$l2Penalty, ncol(regressionIn)-1))
  } else {
    stop("Unknown targetType", opts$targetType)
  }

  outWeightMatrix <- DEEButil::saveSolve(
    XTX,
    crossprod(X, regressionOut))

  return(lst(
    outWeightMatrix,
    timeStep,
    obs,
    valid = TRUE))
}


createFeaturesOneTrajOneTime <- function(traj, row, timeStep, timeStepAsInput, pastSteps, skip, polyDeg = NULL) {
  nRowsRequired <- 1 + pastSteps*(skip+1)
  if (row > nRowsRequired) {
    traj <- traj[(row-nRowsRequired+1):row, ]
  } else {
    traj <- traj[1:row, ]
  }
  if (timeStepAsInput) {
    timeSteps <- c(diff(traj$time), timeStep)
    state <- cbind(timeSteps, traj$state)
  } else {
    state <- traj$state
  }
  linFeatures <- createLinFeaturesOneTrajLastTime(state, pastSteps, skip)
  if (hasValue(polyDeg)) {
    features <- createPolyFeaturesOneTrajOne(linFeatures, polyDeg)
  } else {
    features <- linFeatures
  }
  sel <- is.na(features)
  if (any(sel)) {
    features[sel] <- 0
    cat("Replace NA features by 0.\n")
  }
  return(features)
}


predictLinear <- function(parms, opts, startState, len) {

  opts <- asOpts(opts, c("Linear", "Propagator", "HyperParms"))

  nDims <- length(startState)
  outStates <- matrix(NA_real_, nrow = len+1, ncol = nDims)
  outStates[1, ] <- startState

  if (!parms$valid) return(outStates)

  # Need pastSteps*(skip-1) additional states before startState to start prediction correctly.
  nRowsRequired <- 1 + opts$pastSteps*(opts$skip+1)
  trajPrevious <- NULL
  iStart <- DEEButil::whichMinDist(parms$obs$state, startState)
  if (sum((parms$obs$state[iStart,] - startState)^2) < sqrt(.Machine$double.eps)) {
    trajId <- parms$obs$trajId[iStart]
    traj <- parms$obs[1:iStart, ]
    traj <- traj[traj$trajId == trajId, ]
    trajPrevious <- traj[pmax(1, (nrow(traj)-nRowsRequired+1)):nrow(traj), ]
    cat("Found startState in training data. Use it to initialize features.\n")
    features <- createFeaturesOneTrajOneTime(trajPrevious, nrow(trajPrevious), parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, opts$polyDeg)
  } else {
    cat("Did not find startState in training data. Use startState to initialize features.\n")
    features <- createFeaturesOneTrajOneTime(makeTrajs(time=0, state=outStates[1, , drop=FALSE]), 1, parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, opts$polyDeg)
  }

  prevState <- startState
  for (i in seq_len(len)) {
    prediction <- as.vector(features %*% parms$outWeightMatrix)
    newState <- getPropagatorNextState(prevState, parms$timeStep, prediction, opts$targetType)
    outStates[i+1,] <- newState
    if (any(!is.finite(newState))) break
    trajPrevious$state <- rbind(trajPrevious$state[-1,], newState)
    trajPrevious$time <- c(trajPrevious$time[-1], last(trajPrevious$time)+parms$timeStep) # TODO: time might be strange: have absolute vs need diff time
    features <- createFeaturesOneTrajOneTime(trajPrevious, nrow(trajPrevious), parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, opts$polyDeg)
    prevState <- newState
  }

  return(outStates)
}
