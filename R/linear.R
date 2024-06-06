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


createLinear <- function(obs, timeStepAsInput, pastSteps, skip, polyDeg, l2Penalty) {

  featureSeries <- mapTrajs2Trajs(obs, \(traj) {

    if (timeStepAsInput) {
      timeSteps <- c(diff(traj$time), NA)
    } else {
      timeSteps <- NULL
    }

    featuresLin <- cbind(timeSteps, traj$state)

    for (j in seq_len(pastSteps)) {
      rowIdxs <- seq_len(nrow(traj$state))
      rowIdxs <- rowIdxs - j*(skip+1)
      nEmpty <- sum(rowIdxs <= 0)
      rowIdxs <- rowIdxs[rowIdxs>=1]
      if (timeStepAsInput) {
        featuresLin <- cbind(
          featuresLin,
          c(rep(NA_real_, nEmpty),
            timeSteps[rowIdxs]))
      }
      featuresLin <- cbind(
        featuresLin,
        rbind(
          matrix(NA_real_, nrow = nEmpty, ncol = ncol(traj$state)),
          traj$state[rowIdxs,]))
    }

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

  regressionOut <- do.call(rbind, applyTrajId(obs, \(traj) traj$state[-1,, drop=FALSE]))
  regressionIn <- do.call(
    rbind,
    applyTrajId(featureSeries, \(traj) traj$state[-nrow(traj),, drop=FALSE]))

  naRows <- rowSums(is.na(regressionIn) > 0)
  regressionIn <- regressionIn[!naRows, ]
  regressionOut <- regressionOut[!naRows, ]

  X <- regressionIn
  XTX <- crossprod(X)
  diag(XTX) <- diag(XTX) + c(0, rep(l2Penalty, ncol(regressionIn)))

  outWeightMatrix <- DEEButil::saveSolve(
    XTX,
    crossprod(X, regressionOut))

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE) # mean timeStep

  return(lst(
    outWeightMatrix,
    timeStep,
    timeStepAsInput))
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

  outStates <- matrix(NA_real_, nrow = len+1, ncol = ncol(linear$outWeightMatrix))
  outStates[1, ] <- startState

  # TODO
  stop("Not implemented yet")

  outTrajs <- makeTrajs(
    time = time,
    state = outStates)

  return(outTrajs)
}


predictLinearDeriv <- function(linear, states, derivOrder) {
  t(apply(states, 1, \(s) {
    predictedStates <- predictLinear(linear, s, len = derivOrder)$state
    polyInterpCoeffs <- polynomialInterpolation(esn$timeStep * 0:derivOrder, predictedStates)
    polyInterpCoeffs[2,] # derivative at 0 of polynomial is linear coefficient (second coeff)
  }))
}
