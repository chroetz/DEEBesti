kernConst <- function(delta) {
  rep(1, length(delta))
}

kernRect <- function(delta) {
  (delta > -1 & delta < 1)
}

kernTri <- function(delta) {
  (delta > -1 & delta < 1) * (1-abs(delta))
}

kernParab <- function(delta) {
  (delta > -1 & delta < 1) * (1-delta^2)
}

getKernel <- function(name) {
  switch(
    name,
    Normal = stats::dnorm,
    Const = kernConst,
    Rect = kernRect,
    Tri = kernTri,
    Parab = kernParab,
    stop("Unknown kernel name: ", name)
  )
}

