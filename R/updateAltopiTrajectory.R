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

  if (hasTrajId(trajs)) {
    trajIds <- getTrajIds(trajs)
    hasObs <- lapply(trajIds, \(trajId) {
      trj <- getTrajsWithId(trajs, trajId)
      ob <- getTrajsWithId(obs, trajId)
      timeDist <- outer(trj$time, ob$time, \(x, y) abs(x - y))
      closest <- apply(timeDist, 2, which.min)
      1:length(trj$time) %in% closest # TODO: only works if there is at most one obs per trajs
    }) |> unlist()
  } else {
    timeDist <- outer(trajs$time, obs$time, \(x, y) abs(x - y))
    closest <- apply(timeDist, 2, which.min)
    1:length(trajs$time) %in% closest
  }
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


updateAltopiTrajMultistep <- function(z, obs, gamma, coeff) {

  with(z, {

    A1A1 <- Matrix::bandSparse(
      steps*d, steps*d, c(0, 2),
      diagonals = list(
        c(rep(1, d), rep(2, (steps-2)*d), rep(1, d)),
        rep(-1, (steps-1)*d)),
      symmetric = TRUE) / step_size^2

    if (length(coeff) < 2) {
      b1 <- as.vector(t(a))
    } else {
      coeff <- coeff / sum(coeff)
      W <- Matrix::bandSparse(
        (steps-1)*d, (steps-1)*d,
        -(seq_along(coeff)-1)*2,
        diagonals = c(
          list(c(rep(1, d*(length(coeff)-1)), rep(coeff[1], d*(steps-length(coeff))))),
          lapply(
            2:length(coeff),
            \(j) c(rep(0, (length(coeff)-j) * d), rep(coeff[j], d*(steps-length(coeff))))
          )
        )
      )
      b1 <- W %*% as.vector(t(a))
    }

    A1b1 <- c(
      -b1[1:d],
      b1[1:((steps-2)*d)] - b1[(1+d):((steps-1)*d)],
      b1[((steps-2)*d+1):((steps-1)*d)]) / step_size

    A2A2 <- Matrix::bandSparse(
      steps*d, steps*d, 0,
      diagonals = list(rep(as.numeric(has_obs), each=d)),
      symmetric = TRUE)

    b2 <- double(steps*d)
    b2[rep(has_obs, each=d)] <- as.vector(t(obs$u))

    A <- gamma/steps * A1A1 + (1-gamma)/n * A2A2
    b <- gamma/steps * A1b1 + (1-gamma)/n * b2

    new_trajectory <<- Matrix::solve(A, b)
    return(t(matrix(new_trajectory, nrow = d)))
  })
}
