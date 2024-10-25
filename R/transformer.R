createTransformer <- function(obs, opts) {

  opts <- asOpts(opts, c("Transformer", "Propagator", "HyperParms"))

  tensorflow::set_random_seed(opts$seed)

  contextLen <- opts$contextLen
  stateDim <- ncol(obs$state)
  outState <- obs$state

  timeStep <- getTimeStepTrajs(obs, requireConst=FALSE)

  if (opts$timeStepAsInput) {
    inState <- cbind(outState, c(diff(obs$time), timeStep))
    inDim <- stateDim + 1
  } else {
    inState <- outState
    inDim <- stateDim
  }
  featureDim <- getFeatureDim(inDim, opts)

  posVec <- transformerGetPositionInfo(opts)

  trainSetting <- expand.grid(
    iStart = seq(1, nrow(obs) - contextLen, by = opts$slidingWindowStep))

  xTrain <- sapply(
    seq_len(nrow(trainSetting)),
    \(i) {
      s <- trainSetting[i,,drop=FALSE]
      x <- inState[s$iStart:(s$iStart+contextLen-1),]
      transformerAddPositionInfo(x, posVec)
    })
  dim(xTrain) <- c(contextLen, featureDim, length(xTrain) / (contextLen * featureDim))
  xTrain <- aperm(xTrain, c(3, 1, 2))
  yTrain <- sapply(
    seq_len(nrow(trainSetting)),
    \(i) {
      s <- trainSetting[i,,drop=FALSE]
      x <- outState[s$iStart+contextLen,]
      x
    })
  yTrain <- t(yTrain)

  parms <- trainPropagatorTransformer(xTrain, yTrain, opts)

  return(c(
    parms,
    lst(timeStep, obs)))
}


trainPropagatorTransformer <- function(x, y, opts) {

  model <- buildTransformerModel(
    contextLen = opts$contextLen,
    featureDim = dim(x)[3],
    stateDim = dim(y)[2],
    headSize = opts$headSize,
    nHeads = opts$nHeads,
    nBlocks = opts$nBlocks,
    mlpUnits = opts$mlpUnits,
    mlpDropout = opts$mlpDropout,
    dropout = opts$dropout
  )

  model %>% keras::compile(
    loss = "mse",
    optimizer = keras::optimizer_adam(learning_rate = opts$learningRate)
  )

  callbacks <- list(
    keras::callback_early_stopping(patience = 40, restore_best_weights = TRUE))

  history <-
    model %>%
    keras::fit(
      verbose = 2,
      x,
      y,
      batch_size = opts$batchSize,
      epochs = opts$epochs,
      callbacks = callbacks,
      validation_split = opts$validationSplit,
      shuffle = TRUE
    )

  return(lst(model, history))
}


predictTransformer <- function(parms, opts, startState, len) {

  opts <- asOpts(opts, c("Transformer", "Propagator", "HyperParms"))

  stateDim <- ncol(startState)
  contextLen <- opts$contextLen
  inDim <- stateDim + ifelse(opts$timeStepAsInput, 1, 0)
  featureDim <- getFeatureDim(inDim, opts)

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
    iStart >= contextLen
  ) {
    if (opts$timeStepAsInput) {
      inState <- cbind(parms$obs$state, c(diff(parms$obs$time), parms$timeStep))
    } else {
      inState <- parms$obs$state
    }
    startContext <- inState[(iStart-contextLen+1):iStart, ]
  } else {
    startContext <- matrix(0, nrow = contextLen, ncol = stateDim)
    startContext[contextLen, ] <- startState
    if (opts$timeStepAsInput) {
      startContext <- cbind(startContext, parms$timeStep)
    }
  }

  posVec <- transformerGetPositionInfo(opts)
  startContext <- transformerAddPositionInfo(startContext, posVec)
  x <- startContext

  dim(x) <- c(1, contextLen, featureDim)

  for (i in seq_len(len)) {
    pred <- parms$model %>% predict(x, verbose=2)
    outStates[i+1,] <- pred
    if (any(!is.finite(pred))) break
    x[1, -contextLen, seq_len(inDim)] <- x[1, -1, seq_len(inDim)]
    x[1, contextLen, seq_len(stateDim)] <- pred
    if (opts$timeStepAsInput) {
      x[1, contextLen, inDim] <- parms$timeStep
    }
  }

  return(outStates)
}



buildTransformerModel <- function(
    contextLen,
    featureDim,
    stateDim,
    headSize,
    nHeads = 4,
    nBlocks = 2,
    mlpUnits = c(32),
    mlpDropout = 0.1,
    dropout = 0.1
) {

  inputs <- keras::layer_input(c(contextLen, featureDim))

  x <- inputs
  for (i in 1:nBlocks) {
    x <- x %>%
      transformerEncoder(
        headSize = headSize,
        nHeads = nHeads,
        dropout = dropout
      )
  }

  x <- x %>%
    keras::layer_global_average_pooling_1d(data_format = "channels_first")

  for (dim in mlpUnits) {
    x <- x %>%
      keras::layer_dense(dim, activation = "relu") %>%
      keras::layer_dropout(mlpDropout)
  }

  outputs <- x %>%
    keras::layer_dense(stateDim, activation = "linear")

  model <- keras::keras_model(inputs, outputs)

  return(model)
}


transformerEncoder <- function(
    inputs,
    headSize,
    nHeads,
    dropout = 0
) {

  attentionLayer <- keras::layer_multi_head_attention(
    key_dim = headSize,
    num_heads = nHeads,
    dropout = dropout)

  nFeatures <- dim(inputs) %>% tail(1)

  x <- inputs %>%
    attentionLayer(., .) %>%
    keras::layer_dropout(dropout) %>%
    keras::layer_layer_normalization(epsilon = 1e-6)

  res <- x + inputs

  x <- res %>%
    keras::layer_dense(nFeatures, activation = "relu") |>
    keras::layer_layer_normalization(epsilon = 1e-6)

  return(x + res)
}


transformerGetPositionInfo <- function(opts) {
  switch(
      opts$posEncoding,
      oneHot = diag(opts$contextLen),
      sinusoidal = {
        pos <- seq(0, 1, length.out = opts$contextLen)
        s <- opts$posEncodeFactor*2/opts$posDim*log2(opts$contextLen)
        sapply(
          seq_len(opts$posDim),
          \(i) if(i %% 2 == 0) sin(pi*pos * 2^(s*(i/2-1))) else cos(pi*pos * 2^(s*(i-1)/2)))
        }
    )
}

transformerAddPositionInfo <- function(states, posVec) {
  return(cbind(states, posVec))
}


getFeatureDim <- function(stateDim, opts) {
  stateDim + opts$posDim
}
