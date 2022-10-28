initAltOpt <- function(obs, interSteps) {
  mapTrajs2Trajs(obs, \(trj) {
    n <- getCount(trj)
    steps <- (n-1) * interSteps + 1
    time <- seq(min(trj$time), max(trj$time), length.out = steps)
    trj <- interpolateTrajs(trj, time)
    trj <- setDeriv(trj, "center")
    trj
  })
}


oneAltOptStep <- function(
    trajs, obs, gamma,
    fitDeriv, fitDerivOpts=NULL,
    updateTraj, updateTrajOpts=NULL
  ) {
  newTrajs <- trajs
  newTrajs$deriv <- do.call(fitDeriv, c(list(x=trajs$state, y=trajs$deriv), fitDerivOpts))
  newTrajs$state <- do.call(updateTraj, c(list(traj=newTrajs, obs=obs, gamma=gamma), updateTrajOpts))
  return(newTrajs)
}


fitTrajsAltOpt <- function(obs, hyperParms, memoize = FALSE) {
  hyperParms <- asOpts(hyperParms, c("AltOpt", "HyperParms"))
  if (hyperParms$steps <= 0) return(initAltOpt(obs, hyperParms$interSteps))
  if (memoize) {
    traj <- getFromMemory(hyperParms)
    if (!is.null(traj)) return(traj)
  }
  preHyperParms <- getHyperParmPredecessorAltOpt(hyperParms)
  preTraj <- fitTrajsAltOpt(obs, preHyperParms, memoize)
  nextTraj <- oneAltOptStep(
    preTraj,
    obs,
    gamma = hyperParms$gamma,
    fitDeriv = buildFitter(hyperParms$fitter),
    updateTraj = updateAltOptTraj)
  if (memoize) {
    addToMemory(hyperParms, nextTraj)
  }
  return(nextTraj)
}


getHyperParmPredecessorAltOpt <- function(hyperParms) {
  hyperParms$steps <- hyperParms$steps - 1
  return(hyperParms)
}


updateAltOptTraj <- function(trajs, obs, gamma) {

  stepSize <- trajs$time[2] - trajs$time[1]
  counts <- unname(getCount(trajs))
  countTotal <- sum(counts)
  nTotal <- sum(getCount(obs))

  diagLeft   <- lapply(counts, \(m) c( 0, rep(-1/2, m-2), -1)) |> unlist()
  diagMiddle <- lapply(counts, \(m) c(-1, rep(   0, m-2),  1)) |> unlist()
  diagRight  <- lapply(counts, \(m) c( 1, rep( 1/2, m-2),  0)) |> unlist()
  matDeriv <- Matrix::bandSparse(
    countTotal, countTotal, c(-1,0,1),
    diagonals = list(
      diagLeft[-1],
      diagMiddle,
      diagRight[-countTotal])) / stepSize
  matDerivSymm <- Matrix::crossprod(matDeriv)

  hasObs <- apply2TrajId(
    trajs, obs, simplify = TRUE,
    \(trj, ob) {
      timeDist <- outer(trj$time, ob$time, \(x, y) abs(x - y))
      closest <- apply(timeDist, 2, which.min)
      1:length(trj$time) %in% closest # TODO: only works if there is at most one obs per trajs
  })
  matObsSymm <- Matrix::bandSparse(
    countTotal, countTotal, 0,
    diagonals = list(as.numeric(hasObs)),
    symmetric = TRUE)

  mat <-
    gamma/countTotal * matDerivSymm +
    (1-gamma)/nTotal * matObsSymm

  newState <- sapply(seq_len(getDim(trajs)), \(k) {

    targetDeriv <- as.vector(trajs$deriv[,k])
    targetDerivTransformed <- Matrix::crossprod(matDeriv, targetDeriv)

    tragetObs <- double(countTotal)
    tragetObs[hasObs] <- as.vector(obs$state[,k])

    traget <-
      gamma/countTotal * targetDerivTransformed +
      (1-gamma)/nTotal * tragetObs

    as.matrix(Matrix::solve(mat, traget))
  })

  return(newState)
}



