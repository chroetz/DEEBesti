createRegression <- function(obs, opts) {

  opts <- asOpts(opts, c("Regression", "Propagator", "HyperParms"))

  featureSeries <- createBaseFeatures(obs, opts$timeStepAsInput,  opts$pastSteps, opts$skip)

  regressionOut <- getPropagatorRegressionOut(obs, opts$targetType)
  regressionIn <- do.call(
    rbind,
    applyTrajId(featureSeries, \(traj) traj$state[-nrow(traj),, drop=FALSE]))

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

  nDims <- length(startState)
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
    features <- createFeaturesOneTrajOneTime(trajPrevious, nrow(trajPrevious), parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, polyDeg = NULL)
  } else {
    cat("Did not find startState in training data. Use startState to initialize features.\n")
    features <- createFeaturesOneTrajOneTime(makeTrajs(time=0, state=outStates[1, , drop=FALSE]), 1, parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, polyDeg = NULL)
  }

  prevState <- startState
  for (i in seq_len(len)) {
    prediction <- inferPropagatorRegression(parms, opts, matrix(features, nrow=1)) |> as.vector()
    newState <- getPropagatorNextState(prevState, parms$timeStep, prediction, opts$targetType)
    outStates[i+1,] <- newState
    trajPrevious$state <- rbind(trajPrevious$state[-1,], newState)
    trajPrevious$time <- c(trajPrevious$time[-1], last(trajPrevious$time) + parms$timeStep) # TODO: time might be strange: have absolute vs need diff time
    features <- createFeaturesOneTrajOneTime(trajPrevious, nrow(trajPrevious), parms$timeStep, opts$timeStepAsInput, opts$pastSteps, opts$skip, polyDeg = NULL)
    prevState <- newState
  }

  return(outStates)
}



trainPropagatorRegression <- function(x, y, optsPropagator)  {
  opts <- asOpts(optsPropagator$propagatorRegression, c("PropagatorRegression"))
  name <- ConfigOpts::getClassAt(opts, 2)
  parms <- switch(
    name,
    LocalConst = ,
    LocalLinear = ,
    GaussianProcess = {
      knnFun <- FastKNN::buildKnnFunction(x, opts$neighbors)
      lst(x, y, knnFun)
    },
    NeuralNet = {
      model <- buildNeuralPropagatorModel(opts, ncol(x), ncol(y))
      model %>% keras::compile(
        loss = "mse",
        optimizer = keras::optimizer_adam(learning_rate = opts$learningRate)
      )
      callbacks <- list(
        keras::callback_early_stopping(patience = 100, restore_best_weights = TRUE))
      history <- model %>%
        keras::fit(
          verbose = 2,
          x,
          y,
          batch_size = opts$batchSize,
          epochs = opts$epochs,
          callbacks = callbacks,
          validation_split = opts$validationSplit,
          shuffle = TRUE)
      lst(model, history)
    },
    stop("Unknown name", name)
  )
  return(c(list(dim = ncol(y)), parms))
}


buildNeuralPropagatorModel <- function(opts, inDim, outDim) {
  model <- keras::keras_model_sequential(input_shape = inDim)
  for (layer in opts$layers) {
    model <-
      model %>%
      keras::layer_dense(units = layer, activation = opts$activation)
  }
  model <-
      model %>%
      keras::layer_dense(units = outDim, activation = "linear")
  return(model)
}



inferPropagatorRegression <- function(parms, opts, xout) {
  fit <- matrix(NA_real_, nrow=nrow(xout), ncol=parms$dim)
  for (i in seq_len(nrow(xout))) {
    fit[i, ] <- predictPropagatorRegression(parms, xout[i,, drop=FALSE], opts$propagatorRegression)
  }
  return(fit)
}


predictPropagatorRegression <- function(parms, xout, opts) {
  opts <- asOpts(opts, c("PropagatorRegression"))
  name <- ConfigOpts::getClassAt(opts, 2)
  switch(
    name,
    LocalConst = predictPropagatorRegressionLocalConst(parms, xout, opts),
    LocalLinear = predictPropagatorRegressionLocalLinear(parms, xout, opts),
    GaussianProcess = predictPropagatorRegressionGaussianProcess(parms, xout, opts),
    NeuralNet = predictPropagatorRegressionNeuralNet(parms, xout),
    stop("Unknown name", name)
  )
}


predictPropagatorRegressionGaussianProcess <- function(parms, xout, opts)  {
  knn <- parms$knnFun(xout)
  y <- parms$y[knn$idx, , drop=FALSE]
  x <- parms$x[knn$idx, , drop=FALSE]
  kernelMatrix <- DEEButil::expKernelMatrix(x, opts$bandwidth, opts$regulation)
  kernelVector <- DEEButil::expKernelVectorFromDistSqr(knn$distSqr, opts$bandwidth)
  crossprod(kernelVector, DEEButil::saveSolve(kernelMatrix, y))
}


predictPropagatorRegressionLocalConst <- function(parms, xout, opts) {
  knn <- parms$knnFun(xout)
  y <- parms$y[knn$idx, , drop=FALSE]
  w <- getKernel(opts$kernel)(sqrt(knn$distSqr) / opts$bandwidth)
  yout <- weightedMean(y, w)
  return(yout)
}


predictPropagatorRegressionLocalLinear <- function(parms, xout, opts) {
  knn <- parms$knnFun(xout)
  y <- parms$y[knn$idx, , drop=FALSE]
  x <- parms$x[knn$idx, , drop=FALSE]
  w <- getKernel(opts$kernel)(sqrt(knn$distSqr) / opts$bandwidth)
  if (sum(abs(w)) == 0) w[] <- 1
  X <- cbind(1, x)
  Xw <- X * w
  beta <- DEEButil::saveSolve(crossprod(Xw, X), crossprod(Xw, y))
  crossprod(c(1, xout), beta)
}


predictPropagatorRegressionNeuralNet <- function(parms, xout) {
  parms$model %>% predict(xout, verbose=2)
}

