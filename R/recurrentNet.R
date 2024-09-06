createRecurrentNet <- function(obs, opts) {

  opts <- asOpts(opts, c("RecurrentNet", "Propagator", "HyperParms"))

  tensorflow::set_random_seed(opts$seed) # TODO: use this whenever KERAS is used

  contextLen <- opts$chunkLen
  stateDim <- ncol(obs$state)
  featureDim <- stateDim
  trainState <- obs$state

  trainSetting <- expand.grid(
    iStart = seq(1, nrow(obs) - contextLen, by = opts$slidingWindowStep))

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

  trainModel <- buildRnnModel(opts, dim(x), dim(y), stateful=FALSE)
  trainModel %>% keras::compile(
    loss = "mse",
    optimizer = keras::optimizer_adam(learning_rate = opts$learningRate)
  )
  callbacks <- list(
    keras::callback_early_stopping(patience = 100, restore_best_weights = TRUE))
  history <- trainModel %>%
    keras::fit(
      verbose = 2,
      x,
      y,
      batch_size = opts$batchSize,
      epochs = opts$epochs,
      callbacks = callbacks,
      validation_split = opts$validationSplit,
      shuffle = TRUE)

  return(lst(model = trainModel, history))
}


buildRnnModel <- function(opts, inDim, outDim, stateful) {
  model <- keras::keras_model_sequential(input_shape = inDim[-1], batch_size = if (stateful) 1 else NULL)
  layerFunction <- switch(
    tolower(opts$kind),
    rnn = keras::layer_simple_rnn,
    lstm = keras::layer_lstm,
    gru = keras::layer_gru)

  nLayers <- length(opts$layers)
  for (i in seq_len(nLayers)) {
    layer <- opts$layers[i]
    model <-
      model %>%
      layerFunction(
        units = layer,
        activation = opts$activation,
        stateful = stateful,
        return_sequences = i < nLayers)
  }
  model <-
    model %>%
    keras::layer_dense(units = outDim[-1])
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

  # Decide how to initialize the context
  iStart <- DEEButil::whichMinDist(parms$obs$state, startState)
  if (length(iStart) != 1 || is.na(iStart) || iStart <= 0) {
    warning("Did not find a valid start state. Returning NA.", immediate.=TRUE)
    return(outStates)
  }

  if (
    sum((parms$obs$state[iStart,] - startState)^2) < sqrt(.Machine$double.eps) &&
    iStart >= opts$chunkLen
  ) {
    startContext <- parms$obs$state[(iStart-opts$chunkLen+1):iStart, ]
  } else {
    startContext <- matrix(0, nrow = opts$chunkLen, ncol = stateDim)
    startContext[opts$chunkLen, ] <- startState
  }

  x <- startContext
  dim(x) <- c(1, opts$chunkLen, stateDim)

  for (i in seq_len(len)) {
    pred <- parms$model %>% predict(x, verbose=2)
    outStates[i+1,] <- pred
    if (any(!is.finite(pred))) break
    x[1, -opts$chunkLen, seq_len(stateDim)] <- x[1, -1, seq_len(stateDim)]
    x[1, opts$chunkLen, seq_len(stateDim)] <- pred
  }

  return(outStates)
}
