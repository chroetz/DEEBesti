createRegression <- function(obs, opts) {

  opts <- asOpts(opts, c("Regression", "Propagator", "HyperParms"))

  featureSeries <- createBaseFeatures(obs, opts$timeStepAsInput,  opts$pastSteps, opts$skip)

  regressionOut <- getPropagatorRegressionOut(obs, opts$targetType)
  regressionIn <- do.call(
    rbind,
    applyTrajId(featureSeries$featuresLinTrajs, \(traj) traj$state[-nrow(traj),, drop=FALSE]))

  naRows <- rowSums(is.na(regressionIn) > 0)
  regressionIn <- regressionIn[!naRows, ]
  regressionOut <- regressionOut[!naRows, ]

  parms <- trainPropagatorRegression(regressionIn, regressionOut, opts)

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE)

  return(c(
    parms,
    lst(timeStep, obs)))
}


predictRegression <- function(parms, opts, startState, len) {

  opts <- asOpts(opts, c("Regression", "Propagator", "HyperParms"))

  nDims <- ncol(parms$y)
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
    features <- createFeaturesOne(trajPrevious, nrow(trajPrevious), parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, polyDeg = NULL)
  } else {
    cat("Did not find startState in training data. Use startState to initialize features.\n")
    features <- createFeaturesOne(makeTrajs(time=0, state=outStates[1, , drop=FALSE]), 1, parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, polyDeg = NULL)
  }

  prevState <- startState
  for (i in seq_len(len)) {
    prediction <- inferPropagatorRegression(parms, opts, matrix(features, nrow=1)) |> as.vector()
    newState <- getPropagatorNextState(prevState, parms$timeStep, prediction, opts$targetType)
    outStates[i+1,] <- newState
    trajPrevious$state <- rbind(trajPrevious$state[-1,], newState)
    trajPrevious$time <- c(trajPrevious$time[-1], last(trajPrevious$time) + parms$timeStep) # TODO: time might be strange: have absolute vs need diff time
    features <- createFeaturesOne(trajPrevious, nrow(trajPrevious), parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, polyDeg = NULL)
    prevState <- newState
  }

  return(outStates)
}



trainPropagatorRegression <- function(x, y, opts)  {
  knnFun <- FastKNN::buildKnnFunction(x, opts$neighbors)
  return(lst(x, y, knnFun))
}


inferPropagatorRegression <- function(parms, opts, xout) {
  fit <- matrix(NA_real_, nrow=nrow(xout), ncol=ncol(parms$y))
  for (i in seq_len(nrow(xout))) {
    xouti <- xout[i,]
    knn <- parms$knnFun(xouti)
    x <- parms$x[knn$idx, , drop=FALSE]
    y <- parms$y[knn$idx, , drop=FALSE]
    fit[i, ] <- predictPropagatorRegression(knn$distSqr, x, y, xouti, opts$propagatorRegression)
  }
  return(fit)
}


predictPropagatorRegression <- function(distSqr, x, y, xout, opts) {
  opts <- asOpts(opts, c("PropagatorRegression"))
  name <- ConfigOpts::getClassAt(opts, 2)
  switch(
    name,
    LocalConst = predictPropagatorRegressionLocalConst(distSqr, y, opts),
    LocalLinear = predictPropagatorRegressionLocalLinear(distSqr, x, y, xout, opts),
    GaussianProcess = predictPropagatorRegressionGaussianProcess(distSqr, x, y, opts),
    stop("Unknown name", name)
  )
}


predictPropagatorRegressionLocalConst <- function(distSqr, x, y, opts) {
  w <- getKernel(opts$kernel)(sqrt(distSqr) / opts$bandwidth)
  yout <- weightedMean(y, w)
  return(yout)
}


predictPropagatorRegressionGaussianProcess <- function(distSqr, x, y, opts)  {
  kernelMatrix <- DEEButil::expKernelMatrix(x, opts$bandwidth, opts$regulation)
  kernelVector <- DEEButil::expKernelVectorFromDistSqr(distSqr, opts$bandwidth)
  crossprod(kernelVector, DEEButil::saveSolve(kernelMatrix, y))
}


predictPropagatorRegressionLocalLinear <- function(distSqr, x, y, xout, opts) {
  w <- getKernel(opts$kernel)(sqrt(distSqr) / opts$bandwidth)
  X <- cbind(1, x)
  Xw <- X * w
  beta <- DEEButil::saveSolve(crossprod(Xw, X), crossprod(Xw, y))
  crossprod(c(1, xout), beta)
}

