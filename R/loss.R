#' @export
loss <- function(z, obs, gamma) {
  with(
    z,
    (1-gamma) * mean((obs$u - u[has_obs,])^2) +
      + gamma * mean(((u[2:steps,] - u[1:(steps-1),])/step_size - a)^2)
  )
}

#' @export
err <- function(z, truth) {
  sum((truth$u[1:sum(z$has_obs),] - z$u[z$has_obs,])^2)
}

#' @export
emprisk <- function(z, obs) {
  sum((obs$u - z$u[z$has_obs,])^2)
}

#' @export
change <- function(z_old, z_new) {
  sum((z_old$u - z_new$u)^2/sum(z_old$u^2))
}
