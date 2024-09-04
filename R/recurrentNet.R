createRecurrentNet <- function(obs, opts) {

  opts <- asOpts(opts, c("RecurrentNet", "Propagator", "HyperParms"))

  set.seed(opts$seed) # TODO: set keras/TF seeds

  browser()

  contextLen <- opts$contextLen
  stateDim <- ncol(obs$state)
  featureDim <- stateDim
  trainState <- obs$state

  trainSetting <- expand.grid(
    iStart = seq_len(nrow(obs) - contextLen))

  xTrain <- sapply(
    seq_len(nrow(trainSetting)),
    \(i) {
      s <- trainSetting[i,,drop=FALSE]
      trainState[s$iStart:(s$iStart+contextLen-1),]
    })
  dim(xTrain) <- c(contextLen, featureDim, length(xTrain) / (contextLen * featureDim))
  xTrain <- aperm(xTrain, c(3, 1, 2))
  yTrain <- sapply(
    seq_len(nrow(trainSetting)),
    \(i) {
      s <- trainSetting[i,,drop=FALSE]
      x <- trainState[s$iStart+contextLen,]
      x
    })
  yTrain <- t(yTrain)

  parms <- trainPropagatorRecurrentNet(xTrain, yTrain, opts)

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE)

  return(c(
    parms,
    lst(timeStep, obs)))
}


trainPropagatorRecurrentNet <- function(x, y, opts) {

  # TODO
  name <- ConfigOpts::getClassAt(opts, 4)
  parms <- switch(
    name,
    Lstm = ,
    Rnn = {
      model <- buildRnnModel(opts, dim(x), dim(y))
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


buildRnnModel <- function(opts, inDim, outDim) {
  model <- keras::keras_model_sequential(input_shape = inDim[-1])
  nLayers <- length(opts$layers)
  for (i in seq_len(nLayers)) {
    layer <- opts$layers[i]
    model <-
      model %>%
      keras::layer_simple_rnn(units = layer, activation = opts$activation, return_sequences = i < nLayers)
  }
  return(model)
}


predictRecurrentNet <- function(parms, opts, startState, len) {

  opts <- asOpts(opts, c("RecurrentNet", "Propagator", "HyperParms"))

  stateDim <- length(startState)

  outStates <- matrix(NA_real_, nrow = len+1, ncol = stateDim)
  outStates[1, ] <- startState
  if (any(is.na(startState))) {
    warning("startState has NA. Returning NA.", immediate.=TRUE)
    return(outStates)
  }

  # TODO

  # Decide how to initialize the context
  iStart <- DEEButil::whichMinDist(parms$states, startState)
  if (length(iStart) != 1 || is.na(iStart) || iStart <= 0) {
    warning("Did not find a valid start state. Returning NA.", immediate.=TRUE)
    return(outStates)
  }
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
    # TODO
  }

  return(outStates)
}
