fitTrajsTrajOpt <- function(obs, hyperParms, memoize = FALSE) {
  hyperParms <- asOpts(hyperParms, c("TrajOpt", "FitTrajs"))
  if (hyperParms$steps <= 0) return(initTrajOpt(obs, hyperParms$interSteps))
  if (memoize) {
    traj <- getFromMemory(hyperParms)
    if (!is.null(traj)) return(traj)
  }
  preHyperParms <- getHyperParmPredecessorTrajOpt(hyperParms)
  preTraj <- fitTrajsTrajOpt(obs, preHyperParms, memoize)
  nextTraj <- oneTrajOptStep(
    preTraj,
    obs,
    gamma = hyperParms$gamma,
    neighbors = hyperParms$neighbors,
    derivWeights = hyperParms$derivWeights)
  if (memoize) {
    addToMemory(hyperParms, nextTraj)
  }
  return(nextTraj)
}


initTrajOpt <- function(obs, interSteps) {
  mapTrajs2Trajs(obs, \(trj) {
    n <- getCount(trj)
    steps <- (n-1) * interSteps + 1
    time <- seq(min(trj$time), max(trj$time), length.out = steps)
    trj <- interpolateTrajs(trj, time)
    trj
  })
}



getHyperParmPredecessorTrajOpt <- function(hyperParms) {
  hyperParms$steps <- hyperParms$steps - 1
  return(hyperParms)
}


oneTrajOptStep <- function(trajs, obs, gamma, neighbors, derivWeights) {
  newTrajs <- trajs
  newTrajs$state <- updateTrajOptTraj(trajs, obs, gamma, neighbors, derivWeights)
  return(newTrajs)
}


updateTrajOptTraj <- function(trajs, obs, gamma, neighbors, derivWeights, alpha = 1e-10) {
  # TODO:
  # * multiple trajs
  # * timeSteps <- diff(trajs$time)

  m <- nrow(trajs)
  n <- nrow(obs)
  d <- getDim(trajs)

  knnFun <- FastKNN::buildKnnFunction(trajs$state[-m,], neighbors)
  idxs <- sapply(seq_len(m-1), \(i) {
    knnFun(trajs$state[i,])$idx
  })

  matDeriv <- Matrix::sparseMatrix(
    i = rep(1:((m-1)*neighbors), times=4),
    j = c(
      rep(1:(m-1), each = neighbors),
      rep(2:m, each = neighbors),
      idxs,
      idxs+1),
    x = c(
      rep(-1, (m-1)*neighbors),
      rep(1, (m-1)*neighbors),
      rep(1, (m-1)*neighbors),
      rep(-1, (m-1)*neighbors)),
    dims = c((m-1)*neighbors, m))
  matDerivSymm <- Matrix::crossprod(matDeriv)

  knnFun <- FastKNN::buildKnnFunction(trajs$state[-c(1,m),], neighbors)
  idxs <- sapply(2:(m-1), \(i) {
    knnFun(trajs$state[i,])$idx + 1 # idx refers to trajs$state[-c(1,m),]
  })

  matDeriv2 <- Matrix::sparseMatrix(
    i = rep(1:((m-2)*neighbors), times=6),
    j = c(
      rep(1:(m-2), each = neighbors),
      rep(2:(m-1), each = neighbors),
      rep(3:(m-0), each = neighbors),
      idxs - 1,
      idxs,
      idxs + 1),
    x = c(
      rep(+1, (m-2)*neighbors),
      rep(-2, (m-2)*neighbors),
      rep(+1, (m-2)*neighbors),
      rep(-1, (m-2)*neighbors),
      rep(+2, (m-2)*neighbors),
      rep(-1, (m-2)*neighbors)),
    dims = c((m-2)*neighbors, m))
  matDeriv2Symm <- Matrix::crossprod(matDeriv2)

  hasObs <- apply(
    outer(obs$time, trajs$time, \(a,b) abs(a-b)),
    2,
    \(z) any(z < sqrt(.Machine$double.eps)))
  matObsSymm <- Matrix::bandSparse(
    m, m, 0,
    diagonals = list(as.numeric(hasObs)),
    symmetric = TRUE)

  mat <-
    gamma/m * (derivWeights[1] * matDerivSymm + derivWeights[2] * matDeriv2Symm) +
    (1-gamma)/n * matObsSymm +
    alpha * Matrix::diag(1, nrow = m, ncol = m)

  tragetObs <- sapply(1:d, \(k) {
    y <- double(m)
    y[hasObs] <- as.vector(obs$state[,k])
    y
  })
  traget <-
    (1-gamma)/n * tragetObs +
    alpha * trajs$state


  newState <- as.matrix(Matrix::solve(mat, traget))
  return(newState)
}



