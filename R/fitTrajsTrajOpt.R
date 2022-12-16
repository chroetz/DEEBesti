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
  # * remove derivWeights which are ignored right now

  m <- nrow(trajs)
  n <- nrow(obs)
  d <- getDim(trajs)

  idxs <- t(RANN::nn2(
    trajs$state[-c(1,m),],
    trajs$state[-c(1,m),],
    k = neighbors
  )$nn.idx + 1)
  # write Rcpp function to filter idexes such that a sequence of consecutive indices is reduced to one idx
  # weight elements based on distance and based on uniformity of direction
  # implement a version which does not always include the derivative along the trajectory as one of the two directions

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

  # TODO: make Rcpp implementation of this
  tsDst <- outer(trajs$time, obs$time, \(a,b) abs(a-b))
  iMin <- apply(tsDst, 2, which.min)
  stopifnot(length(unique(iMin)) == length(iMin))
  hasObs <- 1:m %in% iMin

  matObsSymm <- Matrix::bandSparse(
    m, m, 0,
    diagonals = list(as.numeric(hasObs)),
    symmetric = TRUE)

  matOld <- Matrix::bandSparse(
    m, m, 0,
    diagonals = list(rep(alpha, m)),
    symmetric = TRUE)

  mat <-
    gamma/m * matDeriv2Symm +
    (1-gamma)/n * matObsSymm +
    matOld

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



