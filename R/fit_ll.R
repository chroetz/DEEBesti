#' Fit local linear estimator.
#'
#' Take predictor-response pairs and estimate new responses for the same
#' predictors via a local linear estimator. Both, predictors and responses, may
#' be multidimensional.
#'
#' @param x \code{numeric(n)} or \code{matrix(, n, p)}. n observations of
#'   p-dimensional predictors.
#' @param y \code{numeric(n)} or \code{matrix(, n, k)}. n observations of
#'   k-dimensional predictors.
#' @param h \code{numeric(1)}. The bandwidth to use.
#' @param kernel \code{function}. The kernel.
#' @return A numeric n by k matrix. The estimates of response variable at the
#'   locations of the observed predictors.
fit_ll <- function(x, y, h, kernel) {
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  stopifnot(NROW(x) == NROW(y))
  stopifnot(is.numeric(h) && length(h) == 1)
  stopifnot(is.function(kernel))
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  p <- 1 + ncol(x)
  dst <- as.matrix(stats::dist(x))
  Ws <- kernel(dst/h)
  X <- cbind(1, x)
  XTWs <- array(rep(t(X), times = n) * rep(Ws, each = p), dim = c(p, n, n))
  f_hat <- matrix(NA_real_, nrow = n, ncol = ncol(y))
  for (j in 1:n) {
    B <- XTWs[, , j] %*% X
    a <- XTWs[, , j] %*% y
    f_hat[j, ] <- X[j, ] %*% .Internal(La_solve(B, a, .Machine$double.eps))
  }
  return(f_hat)
}
