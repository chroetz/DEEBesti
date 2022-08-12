#' @export
kern_parab <- function(delta) {
  (delta > -1 & delta < 1) * (1-delta^2)
}

#' @export
kern_rect <- function(delta) {
  (delta > -1 & delta < 1)
}

#' @export
kern_tri <- function(delta) {
  (delta > -1 & delta < 1) * (1-abs(delta))
}
