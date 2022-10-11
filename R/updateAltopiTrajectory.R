updateAltopiTraj <- function(trajs, obs, gamma) {

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
