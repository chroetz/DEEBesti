#' @export
update_trajectory <- function(z, obs, gamma) {
  with(z, {

    A1A1 <- Matrix::bandSparse(
      steps*d, steps*d, c(0, 2),
      diagonals = list(
        c(rep(1, d), rep(2, (steps-2)*d), rep(1, d)),
        rep(-1, (steps-1)*d)),
      symmetric = TRUE) / step_size^2

    b1 <- as.vector(t(a))
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


#' @export
update_trajectory_multistep <- function(z, obs, gamma, coeff) {

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
