
parab <- function(delta) {
  (delta > -1 & delta < 1) * (1-delta^2)
}

rect <- function(delta) {
  (delta > -1 & delta < 1)
}

tri <- function(delta) {
  (delta > -1 & delta < 1) * (1-abs(delta))
}
