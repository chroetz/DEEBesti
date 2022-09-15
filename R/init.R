#' @export
init <- function(obs, inter_steps) {
  n <- nrow(obs$state)
  inter_steps <- inter_steps
  steps <- (n-1)*inter_steps + 1
  d <- ncol(obs$state)
  t <- seq(min(obs$time), max(obs$time), length.out = steps)
  step_size <- t[2] - t[1]
  state <- do.call(
    cbind,
    lapply(seq_len(d), \(j) stats::approx(obs$time, obs$state[,j], t)$y))
  deriv <- do.call(
    cbind,
    lapply(seq_len(d), \(j) (state[2:steps,j]-state[1:(steps-1),j]) / step_size))
  has_obs <- 1:steps %% inter_steps == 1
  lst <- as.list(environment())
  return(lst[setdiff(names(lst), "obs")])
}
