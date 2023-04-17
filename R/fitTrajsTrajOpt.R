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
    alpha = hyperParms$alpha)
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


oneTrajOptStep <- function(trajs, obs, gamma, neighbors, alpha) {
  newTrajs <- trajs
  newTrajs$state <- updateTrajOptTraj(trajs, obs, gamma, neighbors, alpha)
  return(newTrajs)
}


updateTrajOptTraj <- function(trajs, obs, gamma, neighbors, alpha) {
  # weight elements based on distance and based on uniformity of direction
  # implement a version which does not always include the derivative along the trajectory as one of the two directions
  PT <- FALSE
  if (PT)cat("Start updateTrajOptTraj()\n")
  pt0 <- proc.time()

  m <- nrow(trajs)
  n <- nrow(obs)
  d <- getDim(trajs)

  isInnerId <- \(x) {
    len <- length(x)
    c(FALSE, x[-1] == x[-len]) & c(x[-1] == x[-len], FALSE)
  }

  isInner <- isInnerId(trajs$trajId)
  idxOffset <- cumsum(!isInner)[isInner]
  innerIdx <- which(isInner)
  nInner <- length(innerIdx)
  innerState <- trajs$state[isInner, ]

  if (PT)cat("NN\n")
  nn <- RANN::nn2(
    innerState,
    innerState,
    k = neighbors + 1
  )

  neighborDists <- t(nn$nn.dists[,-1])
  neighborInnerIdx <- as.vector(nn$nn.idx[,-1])
  neighborIdx <- matrix(
    neighborInnerIdx + idxOffset[neighborInnerIdx],
    nrow = neighbors,
    byrow=TRUE)

  pt1 <- proc.time()
  if (PT)print(pt1-pt0)

  if (PT)cat("sparseMat\n")
  matDeriv2 <- Matrix::sparseMatrix(
    i = rep(1:(nInner*neighbors), times=6),
    j = c(
      rep(innerIdx-1, each = neighbors),
      rep(innerIdx, each = neighbors),
      rep(innerIdx+1, each = neighbors),
      neighborIdx - 1,
      neighborIdx,
      neighborIdx + 1),
    x = c(
      rep(+1, nInner*neighbors)/neighborDists,
      rep(-2, nInner*neighbors)/neighborDists,
      rep(+1, nInner*neighbors)/neighborDists,
      rep(-1, nInner*neighbors)/neighborDists,
      rep(+2, nInner*neighbors)/neighborDists,
      rep(-1, nInner*neighbors)/neighborDists),
    dims = c(nInner*neighbors, m))
  matDeriv2Symm <- Matrix::crossprod(matDeriv2)
  pt2 <- proc.time()
  if (PT)print(pt2-pt1)

  if (PT)cat("time\n")
  subSteps <- round(getTimeStepTrajs(obs)/getTimeStepTrajs(trajs))
  lensObs <- getCount(obs)
  lensTrajs <- getCount(trajs)
  startIdxTrajs <- 1 + c(0, cumsum(lensTrajs[-length(lensTrajs)]))
  iMin <- unlist(lapply(
    seq_along(lensObs),
    \(i) startIdxTrajs[i] + (seq_len(lensObs[i])-1)*subSteps))
  stopifnot(max(abs(trajs$time[iMin] - obs$time)) < sqrt(.Machine$double.eps))
  hasObs <- 1:m %in% iMin
  pt3 <- proc.time()
  if (PT)print(pt3-pt2)

  if (PT)cat("createMats\n")
  matObsSymm <- Matrix::bandSparse(
    m, m, 0,
    diagonals = list(as.numeric(hasObs)),
    symmetric = TRUE)

  matOld <- Matrix::bandSparse(
    m, m, 0,
    diagonals = list(rep(alpha, m)),
    symmetric = TRUE)

  mat <-
    gamma / neighbors * matDeriv2Symm +
    matObsSymm + matOld

  targetObs <- sapply(1:d, \(k) {
    y <- double(m)
    y[hasObs] <- as.vector(obs$state[,k])
    y
  })
  target <-
    targetObs + alpha * trajs$state
  pt4 <- proc.time()
  if (PT)print(pt4-pt3)

  if (PT)cat("solve\n")
  newState <- as.matrix(Matrix::solve(mat, target))
  pt5 <- proc.time()
  if (PT)print(pt5-pt4)

  return(newState)
}



