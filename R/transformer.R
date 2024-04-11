createTransformer <- function(opts, stateDim) {

  posDim <- opts$contextLen
  featureDim <- stateDim + posDim
  posVec <- diag(contextLen)
  dim(posVec) <- c(posDim, contextLen)
  posVec <- t(posVec)

  model <- buildTransformerModel(
    contextLen = opts$contextLen,
    featureDim = featureDim,
    stateDim = stateDim,
    head_size = opts$headSize,
    nHeads = opts$nHeads,
    nBlocks = opts$nBlocks,
    mlpUnits = opts$mlpUnits,
    mlpDropout = opts$mlpDropout,
    dropout = opts$dropout
  )

  model %>% compile(
    loss = "mse",
    optimizer = optimizer_adam(learning_rate = opts$learningRate)
  )

  return(lst(model, posVec))

}

trainTransformer <- function(transformer, obs, opts) {

  trainState <- obs$state

  trainSetting <- expand.grid(
    iStart = seq_len(nrow(obs) - opts$contextLen))

  xTrain <- sapply(
    seq_len(nrow(trainSetting)),
    \(i) {
      s <- trainSetting[i,,drop=FALSE]
      x <- cbind(trainState[s$iStart:(s$iStart+contextLen-1),], transformer$posVec)
      x
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
  shuffledIndex <- sample(seq_len(nrow(yTrain)))
  xTrain <- xTrain[shuffledIndex,,]
  yTrain <- yTrain[shuffledIndex,]

  callbacks <- list(
    keras::callback_early_stopping(patience = 40, restore_best_weights = TRUE))

  history <- transformer$model %>%
    fit(
      xTrain,
      yTrain,
      steps_per_epoch = 100,
      epochs = 500,
      callbacks = callbacks,
      validation_split = 0.1,
      shuffle = TRUE
    )

  return(transformer)
}

predictTransformer <- function(transformer, startState, len = NULL, startTime = 0, timeRange = NULL) {

  stateDim <- ncol(startState)
  posDim <- opts$contextLen
  featureDim <- stateDim + posDim

  # TODO
  firstTest <- matrix(0, nrow=opts$contextLen, ncol=featureDim)
  lastTrainPart <- trainState[(nrow(trainState)-contextLen+1):nrow(trainState),]
  firstTest[, seq_len(stateDim)] <- lastTrainPart
  firstTest[, stateDim + seq_len(posDim)] <- posVec
  dim(firstTest) <- c(1, contextLen, featureDim)
  x <- firstTest
  nTest <- nrow(test)*2
  predStates <- matrix(NA_real_, nrow=nTest, ncol=stateDim)
  for (i in 1:nTest) {
    cat(i, ": ")
    pred <- model %>% predict(x)
    predStates[i,] <- pred
    x[1, -contextLen, seq_len(stateDim)] <- x[1, -1, seq_len(stateDim)]
    x[1, contextLen, seq_len(stateDim)] <- pred
  }

  trainTimeStep <- mean(diff(train$time))

  esti <- DEEBtrajs::makeTrajs(train$time[nrow(train)] + seq_len(nTest)*trainTimeStep, predStates)
  esti <- DEEBtrajs::interpolateTrajs(esti, test$time)

  return(esti)
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
        head_size = headSize,
        num_heads = nHeads,
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
  # Attention and Normalization
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

  # Feed Forward Part
  x <- res %>%
    keras::layer_dense(nFeatures, activation = "relu") |>
    keras::layer_layer_normalization(epsilon = 1e-6)

  # return output + residual
  return(x + res)
}
