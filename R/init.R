#' @export
init <- function(obs, inter_steps) {
  n <- nrow(obs$u)
  inter_steps <- inter_steps
  steps <- (n-1)*inter_steps + 1
  d <- ncol(obs$u)
  t <- seq(min(obs$t), max(obs$t), length.out = steps)
  step_size <- t[2] - t[1]
  u <- do.call(
    cbind,
    lapply(seq_len(d), \(j) stats::approx(obs$t, obs$u[,j], t)$y))
  a <- do.call(
    cbind,
    lapply(seq_len(d), \(j) (u[2:steps,j]-u[1:(steps-1),j]) / step_size))
  has_obs <- 1:steps %% inter_steps == 1
  lst <- as.list(environment())
  return(lst[setdiff(names(lst), "obs")])
}
