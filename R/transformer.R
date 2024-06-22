createTransformer <- function(opts, stateDim) {

  opts <- asOpts(opts, c("Transformer", "Propagator", "HyperParms"))

  featureDim <- getFeatureDim(stateDim, opts)

  model <- buildTransformerModel(
    contextLen = opts$contextLen,
    featureDim = featureDim,
    stateDim = stateDim,
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

  env <- new.env(parent=emptyenv())

  return(lst(model, env))

}


trainTransformer <- function(parms, obs, opts) {

  opts <- asOpts(opts, c("Transformer", "Propagator", "HyperParms"))

  contextLen <- opts$contextLen
  stateDim <- ncol(obs$state)
  featureDim <- getFeatureDim(stateDim, opts)
  trainState <- obs$state

  trainSetting <- expand.grid(
    iStart = seq_len(nrow(obs) - contextLen))

  xTrain <- sapply(
    seq_len(nrow(trainSetting)),
    \(i) {
      s <- trainSetting[i,,drop=FALSE]
      states <- trainState[s$iStart:(s$iStart+contextLen-1),]
      transformerAddPositionInfo(parms, states, opts)
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

  callbacks <- list(
    keras::callback_early_stopping(patience = 40, restore_best_weights = TRUE))

  history <- parms$model %>%
    keras::fit(
      verbose = 2,
      xTrain,
      yTrain,
      batch_size = opts$batchSize,
      epochs = opts$epochs,
      callbacks = callbacks,
      validation_split = opts$validationSplit,
      shuffle = TRUE
    )

  parms$timeStep <- getTimeStepTrajs(obs, requireConst=FALSE) # mean timeStep
  parms$history <- history
  parms$states <- trainState

  return(parms)
}


predictTransformer <- function(parms, opts, startState, len) {

  opts <- asOpts(opts, c("Transformer", "Propagator", "HyperParms"))

  stateDim <- ncol(startState)
  contextLen <- opts$contextLen
  featureDim <- getFeatureDim(stateDim, opts)

  # Decide how to initialize the context
  iStart <- DEEButil::whichMinDist(parms$states, startState)
  if (
    sum((parms$states[iStart,] - startState)^2) < sqrt(.Machine$double.eps) &&
    iStart >= contextLen
  ) {
    startContext <- parms$states[(iStart-contextLen+1):iStart, ]
  } else {
    startContext <- matrix(0, nrow = contextLen, ncol = featureDim)
    startContext[contextLen, ] <- startState
  }

  startContext <- transformerAddPositionInfo(parms, startContext, opts)
  x <- startContext
  dim(x) <- c(1, contextLen, featureDim)
  outStates <- matrix(NA_real_, nrow = len+1, ncol = stateDim)
  outStates[1, ] <- startState

  for (i in seq_len(len)) {
    pred <- parms$model %>% predict(x, verbose=2)
    outStates[i+1,] <- pred
    x[1, -contextLen, seq_len(stateDim)] <- x[1, -1, seq_len(stateDim)]
    x[1, contextLen, seq_len(stateDim)] <- pred
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


transformerAddPositionInfo <- function(transformer, states, opts) {
  if (length(transformer$env$posVec) == 0) {
    transformer$env$posVec <- switch(
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
  return(cbind(states, transformer$env$posVec))
}


getFeatureDim <- function(stateDim, opts) {
  stateDim + opts$posDim
}
