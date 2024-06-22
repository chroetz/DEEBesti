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


createLinFeaturesOne <- function(state, pastSteps, skip) {
  featuresLin <- state
  for (j in seq_len(pastSteps)) {
    rowIdxs <- seq_len(nrow(state))
    rowIdxs <- rowIdxs - j*(skip+1)
    nEmpty <- sum(rowIdxs <= 0)
    rowIdxs <- rowIdxs[rowIdxs>=1]
    featuresLin <- cbind(
      featuresLin,
      rbind(
        matrix(NA_real_, nrow = nEmpty, ncol = ncol(state)),
        state[rowIdxs, , drop=FALSE]))
  }
  return(featuresLin)
}


createBaseFeatures <- function(trajs, timeStepAsInput,  pastSteps, skip) {

  featuresLinTrajs <- mapTrajs2Trajs(trajs, \(traj) {

    featuresLin <- createLinFeaturesOne(trajs$state, pastSteps, skip)

    out <- makeTrajs(
      time = traj$time,
      state = featuresLin)

    return(out)

  })

  if (!isFALSE(timeStepAsInput) && length(timeStepAsInput) > 0) {

     featuresTimeTrajs <- mapTrajs2Trajs(trajs, \(traj) {

      if (is.numeric(timeStepAsInput)) {
        timeSteps <- rep(timeStepAsInput, nrow(traj$state))
      } else if (is.logical(timeStepAsInput) && timeStepAsInput) {
        timeSteps <- c(diff(traj$time), NA)
      } else {
        stop("timeStepAsInput must be logical or numeric")
      }

      featuresTime <- matrix(timeSteps, ncol=1)

      for (j in seq_len(pastSteps)) {
        rowIdxs <- seq_len(nrow(traj$state))
        rowIdxs <- rowIdxs - j*(skip+1)
        nEmpty <- sum(rowIdxs <= 0)
        rowIdxs <- rowIdxs[rowIdxs>=1]
        featuresTime <- cbind(
          featuresTime,
          c(rep(NA_real_, nEmpty),
            timeSteps[rowIdxs]))
      }

      out <- makeTrajs(
        time = traj$time,
        state = featuresTime)

      return(out)

    })

  } else {

    featuresTimeTrajs <- NULL

  }

  return(lst(featuresLinTrajs, featuresTimeTrajs))
}


createPolyFeaturesOne <- function(featuresLin, polyDeg) {
  nCols <- ncol(featuresLin)
  features <- matrix(1, nrow = nrow(featuresLin), ncol=1)
  for (deg in seq_len(polyDeg)) {
    indices <- do.call(expand.grid, replicate(deg, seq_len(nCols), simplify=FALSE))
    indices <- as.matrix(indices)
    notOrdered <- rowSums(indices[, -ncol(indices), drop=FALSE] > indices[, -1, drop=FALSE]) > 0
    indices <- indices[!notOrdered, , drop=FALSE]
    values <- unlist(apply(indices, 1, \(idx) rowIdxProd(featuresLin, idx), simplify=FALSE))
    features <- cbind(features, matrix(values, nrow = nrow(features)))
  }
  return(features)
}


createPolyFeatures <- function(baseFeatures, polyDeg) {

  if (hasValue(baseFeatures$featuresTimeTrajs)) {
    trajs <- baseFeatures$featuresLinTrajs
  } else {
    trajs <- makeTrajs(
      baseFeatures$featuresLinTrajs$time,
      cbind(baseFeatures$featuresLinTrajs$state, baseFeatures$featuresTimeTrajs$state))
  }

  mapTrajs2Trajs(trajs, \(traj) {

    features <- createPolyFeaturesOne(traj$state, polyDeg)

    out <- makeTrajs(
      time = traj$time,
      state = features)

    return(out)

  })
}


createLinear <- function(obs, opts) {

  opts <- asOpts(opts, c("Linear", "Propagator", "HyperParms"))

  baseFeatures <- createBaseFeatures(obs, opts$timeStepAsInput,  opts$pastSteps, opts$skip)
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

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE)

  return(lst(
    outWeightMatrix,
    timeStep,
    obs))
}


createFeaturesOne <- function(traj, row, timeStep, timeStepAsInput, pastSteps, skip, polyDeg = NULL) {
  nRowsRequired <- 1 + pastSteps*(skip+1)
  if (row > nRowsRequired) {
    traj <- traj[(row-nRowsRequired+1):row, ]
  } else {
    traj <- traj[1:row, ]
  }
  timeSteps <- c(diff(traj$time), timeStep)
  linFeatures <- createLinFeaturesOne(traj$state, pastSteps, skip)
  if (timeStepAsInput) {
    linFeatures <- cbind(linFeatures, timeSteps)
  }
  if (hasValue(polyDeg)) {
    features <- createPolyFeaturesOne(linFeatures, polyDeg)
  } else {
    features <- linFeatures
  }
  featuresOne <- features[nrow(features), ]
  sel <- is.na(featuresOne)
  if (any(sel)) {
    featuresOne[sel] <- 0
    cat("Replace NA features by 0.\n")
  }
  return(featuresOne)
}


predictLinear <- function(parms, opts, startState, len) {

  opts <- asOpts(opts, c("Linear", "Propagator", "HyperParms"))

  nDims <- ncol(parms$outWeightMatrix)
  outStates <- matrix(NA_real_, nrow = len+1, ncol = nDims)
  outStates[1, ] <- startState

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
    features <- createFeaturesOne(trajPrevious, nrow(trajPrevious), parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, opts$polyDeg)
  } else {
    cat("Did not find startState in training data. Use startState to initialize features.\n")
    features <- createFeaturesOne(makeTrajs(time=0, state=outStates[1, , drop=FALSE]), 1, parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, opts$polyDeg)
  }

  prevState <- startState
  for (i in seq_len(len)) {
    prediction <- crossprod(parms$outWeightMatrix, features) |> as.vector()
    newState <- getPropagatorNextState(prevState, parms$timeStep, prediction, opts$targetType)
    outStates[i+1,] <- newState
    trajPrevious$state <- rbind(trajPrevious$state[-1,], newState)
    trajPrevious$time <- c(trajPrevious$time[-1], last(trajPrevious$time)+parms$timeStep) # TODO: time might be strange: have absolute vs need diff time
    features <- createFeaturesOne(trajPrevious, nrow(trajPrevious), parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, opts$polyDeg)
    prevState <- newState
  }

  return(outStates)
}
