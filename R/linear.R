rowIdxProd <- function(mat, idxs) {
  n <- nrow(mat)
  if (n == 0) return(double(0))
  if (length(idxs) == 0) return(rep(1, n))
  if (length(idxs) == 1) return(mat[,idxs])
  y <- mat[,idxs[[1]]]
  for (i in idxs[-1]) {
    y <- y * mat[,idxs[[i]]]
  }
  return(y)
}


createBaseFeatures <- function(trajs, timeStepAsInput,  pastSteps, skip) {

  featuresLinTrajs <- mapTrajs2Trajs(trajs, \(traj) {

    featuresLin <- traj$state

    for (j in seq_len(pastSteps)) {
      rowIdxs <- seq_len(nrow(traj$state))
      rowIdxs <- rowIdxs - j*(skip+1)
      nEmpty <- sum(rowIdxs <= 0)
      rowIdxs <- rowIdxs[rowIdxs>=1]
      featuresLin <- cbind(
        featuresLin,
        rbind(
          matrix(NA_real_, nrow = nEmpty, ncol = ncol(traj$state)),
          traj$state[rowIdxs,]))
    }

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



createPolyFeatures <- function(baseFeatures, polyDeg) {

  if (is.null(baseFeatures$featuresTimeTrajs)) {
    trajs <- baseFeatures$featuresLinTrajs
  } else {
    trajs <- makeTrajs(
      baseFeatures$featuresLinTrajs$time,
      cbind(baseFeatures$featuresLinTrajs$state, baseFeatures$featuresTimeTrajs$state))
  }

  mapTrajs2Trajs(trajs, \(traj) {

    featuresLin <- traj$state

    nCols <- ncol(featuresLin)
    features <- matrix(1, nrow = nrow(featuresLin), ncol=1)
    for (deg in polyDeg) {
      indices <- do.call(expand.grid, replicate(deg, seq_len(nCols), simplify=FALSE))
      indices <- as.matrix(indices)
      notOrdered <- rowSums(indices[, -ncol(indices), drop=FALSE] > indices[, -1, drop=FALSE]) > 0
      indices <- indices[!notOrdered, , drop=FALSE]
      features <- cbind(features, apply(indices, 1, \(idx) rowIdxProd(featuresLin, idx)))
    }

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
  diag(XTX) <- diag(XTX) + c(0, rep(l2Penalty, ncol(regressionIn)))

  outWeightMatrix <- DEEButil::saveSolve(
    XTX,
    crossprod(X, regressionOut))

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE) # mean timeStep
  if (is.null(baseFeatures$featuresTimeTrajs)) {
    featureSeriesPredict <- featureSeries
  } else {
    baseFeatures$featuresTimeTrajs$state[is.na(baseFeatures$featuresTimeTrajs$state)] <- timeStep
    featureSeriesPredict <- createPolyFeatures(baseFeatures, polyDeg)
  }

  return(lst(
    timeStepAsInput, pastSteps, skip, polyDeg, l2Penalty,
    outWeightMatrix,
    timeStep,
    features = featureSeriesPredict$state,
    states = obs$state))
}


createFeaturesOne <- function(states, timeStep, timeStepAsInput, pastSteps, skip, polyDeg) {
  validRows <- which(rowSums(is.na(states)) == 0)
  states <- states[seq_len(validRows[length(states)]), ]
  if (nrow(states) > 1+ pastSteps*skip) {
    states <- states[(nrow(states)-pastSteps*skip):nrow(states), ]
  }
  n <- nrow(states)
  trajs <- makeTrajs(
    time = (0:(n-1)) * timeStep,
    state = states)
  baseFeatures <- createBaseFeatures(trajs, timeStepAsInput,  pastSteps, skip)
  featureSeries <- createPolyFeatures(baseFeatures, polyDeg)
  return(featureSeries$state[n, ])
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

  # Need pastSteps*skip additional states before startState to start prediction correctly.
  iStart <- DEEButil::whichMinDist(linear$states, startState)
  if (sum((linear$states[iStart,] - startState)^2) < sqrt(.Machine$double.eps)) {
    cat("Found startState in training data. Use it to initialize features.\n")
    startFeatures <- linear$features[iStart, ]
    features <- startFeatures
  } else {
    cat("Did not find startState in training data. Use startState to initialize features.\n")
    features <- createFeaturesOne(outStates, linear$timeStep, linear$timeStepAsInput, linear$pastSteps, linear$skip, linear$polyDeg)
  }
  sel <- is.na(features)
  if (any(sel)) {
    features[sel] <- 0
    cat("Replace NA start features by 0.\n")
  }

  for (i in seq_len(len)) {
    x <- crossprod(linear$outWeightMatrix, features)
    outStates[i+1,] <- x
    features <- createFeaturesOne(outStates, linear$timeStep, linear$timeStepAsInput, linear$pastSteps, linear$skip, linear$polyDeg)
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
