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


createLinear <- function(obs, timeStepAsInput, pastSteps, skip, polyDeg, l2Penalty) {

  baseFeatures <- createBaseFeatures(obs, timeStepAsInput,  pastSteps, skip)
  featureSeries <- createPolyFeatures(baseFeatures, polyDeg)

  regressionOut <- do.call(rbind, applyTrajId(obs, \(traj) traj$state[-1,, drop=FALSE]))
  regressionIn <- do.call(
    rbind,
    applyTrajId(featureSeries, \(traj) traj$state[-nrow(traj),, drop=FALSE]))

  naRows <- rowSums(is.na(regressionIn) > 0)
  regressionIn <- regressionIn[!naRows, ]
  regressionOut <- regressionOut[!naRows, ]

  # TODO: elastic net: glmnet::glmnet

  X <- regressionIn
  XTX <- crossprod(X)
  diag(XTX) <- diag(XTX) + c(0, rep(l2Penalty, ncol(regressionIn)-1))

  outWeightMatrix <- DEEButil::saveSolve(
    XTX,
    crossprod(X, regressionOut))

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE) # mean timeStep
  if (hasValue(baseFeatures$featuresTimeTrajs)) {
    featureSeriesPredict <- featureSeries
  } else {
    baseFeatures$featuresTimeTrajs$state[is.na(baseFeatures$featuresTimeTrajs$state)] <- timeStep
    featureSeriesPredict <- createPolyFeatures(baseFeatures, polyDeg)
  }

  return(lst(
    timeStepAsInput, pastSteps, skip, polyDeg, l2Penalty,
    outWeightMatrix,
    timeStep,
    obs))
}


createFeaturesOne <- function(traj, row, timeStep, timeStepAsInput, pastSteps, skip, polyDeg) {
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
  features <- createPolyFeaturesOne(linFeatures, polyDeg)
  featuresOne <- features[nrow(features), ]
  sel <- is.na(featuresOne)
  if (any(sel)) {
    featuresOne[sel] <- 0
    cat("Replace NA features by 0.\n")
  }
  return(featuresOne)
}


predictLinear <- function(linear, startState, len = NULL, startTime = 0, timeRange = NULL) {

  if (is.null(timeRange)) {
    stopifnot(length(len) == 1, len >= 0)
    time <- startTime + (0:len)*linear$timeStep
  } else {
    stopifnot(length(timeRange) == 2)
    time <- seq(timeRange[1], timeRange[2], by = linear$timeStep)
    if (time[length(time)] < timeRange[2]) {
      time <- c(time, time[length(time)] + linear$timeStep)
    }
    len <- length(time) - 1
  }

  nDims <- ncol(linear$outWeightMatrix)
  outStates <- matrix(NA_real_, nrow = len+1, ncol = nDims)
  outStates[1, ] <- startState

  # Need pastSteps*(skip-1) additional states before startState to start prediction correctly.
  nRowsRequired <- 1 + linear$pastSteps*(linear$skip+1)
  trajPrevious <- NULL
  iStart <- DEEButil::whichMinDist(linear$obs$state, startState)
  if (sum((linear$obs$state[iStart,] - startState)^2) < sqrt(.Machine$double.eps)) {
    trajId <- linear$obs$trajId[iStart]
    traj <- linear$obs[1:iStart, ]
    traj <- traj[traj$trajId == trajId, ]
    trajPrevious <- traj[pmax(1, (nrow(traj)-nRowsRequired+1)):nrow(traj), ]
    cat("Found startState in training data. Use it to initialize features.\n")
    features <- createFeaturesOne(trajPrevious, nrow(trajPrevious), linear$timeStep, linear$timeStepAsInput, linear$pastSteps, linear$skip, linear$polyDeg)
  } else {
    cat("Did not find startState in training data. Use startState to initialize features.\n")
    features <- createFeaturesOne(makeTrajs(time=0, state=outStates[1, , drop=FALSE]), 1, linear$timeStep, linear$timeStepAsInput, linear$pastSteps, linear$skip, linear$polyDeg)
  }

  for (i in seq_len(len)) {
    x <- crossprod(linear$outWeightMatrix, features) |> as.vector()
    outStates[i+1,] <- x
    trajPrevious$state <- rbind(trajPrevious$state[-1,], x)
    trajPrevious$time <- c(trajPrevious$time[-1], time[i+1]) # TODO: time might be strange: have absolute vs need diff time
    features <- createFeaturesOne(trajPrevious, nrow(trajPrevious), linear$timeStep, linear$timeStepAsInput, linear$pastSteps, linear$skip, linear$polyDeg)
  }

  outTrajs <- makeTrajs(
    time = time,
    state = outStates)

  return(outTrajs)
}


predictLinearDeriv <- function(linear, states, derivOrder) {
  t(apply(states, 1, \(s) {
    predictedStates <- predictLinear(linear, s, len = derivOrder)$state
    polyInterpCoeffs <- polynomialInterpolation(linear$timeStep * 0:derivOrder, predictedStates)
    polyInterpCoeffs[2,] # derivative at 0 of polynomial is linear coefficient (second coeff)
  }))
}
