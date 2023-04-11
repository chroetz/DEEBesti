derivFunGlobalLm <- function(u, parms) {
  d <- length(u)
  z <- double(d)
  for (j in seq_len(d)) {
    fea <- parms$lmFuns$vector$features(u, j)
    z[j] <- sum(parms$coef[[j]] * fea)
  }
  du <- parms$lmFuns$vector$invTransform(u, z)
  return(du)
}


# TODO: document
buildLmFuns <- function(opts) {

  opts <- asOpts(opts, c("GlobalLm", "DerivFun"))
  out <- list(vector = list(), matrix = list())
  # TODO: remove code duplications
  # note: this ignores names of the list entries

  feaMat <-
    opts$features |>
    vapply(
      \(s) gsub("\\bx(\\d+)\\b", "x[,\\1]", s) |> paste0(collapse=", "),
      character(1)
    ) |>
    paste0("cbind(", ...=_, ")") |>
    str2expression()
  if (length(feaMat) == 1) {
    out$matrix$features <- function(x, j) eval(feaMat)
  } else {
    out$matrix$features <- function(x, j) eval(feaMat[[j]])
  }

  feaVec <-
    opts$features |>
    vapply(
      \(s) gsub("\\bx(\\d+)\\b", "x[\\1]", s) |> paste0(collapse=", "),
      character(1)
    ) |>
    paste0("c(", ...=_, ")") |>
    str2expression()
  if (length(feaVec) == 1) {
    out$vector$features <- function(x, j) eval(feaVec)
  } else {
    out$vector$features <- function(x, j) eval(feaVec[[j]])
  }

  if (length(opts$transform) == 0) {
    out$matrix$transform <- function(x, y) y
    out$vector$transform <- function(x, y) y
  } else {
    trnsMat <-
      opts$transform |>
      as.character() |>
      gsub("\\bx(\\d+)\\b", "x[,\\1]", x=_) |>
      gsub("\\by(\\d+)\\b", "y[,\\1]", x=_) |>
      paste0(collapse=", ") |>
      paste0("cbind(", ...=_, ")") |>
      str2expression()
    out$matrix$transform <- function(x, y) eval(trnsMat)
    trnsVec <-
      opts$transform |>
      as.character() |>
      gsub("\\bx(\\d+)\\b", "x[\\1]", x=_) |>
      gsub("\\by(\\d+)\\b", "y[\\1]", x=_) |>
      paste0(collapse=", ") |>
      paste0("c(", ...=_, ")") |>
      str2expression()
    out$vector$transform <- function(x, y) eval(trnsVec)
  }

  if (length(opts$invTransform) == 0) {
    out$matrix$invTransform <- function(x, z) z
    out$vector$invTransform <- function(x, z) z
  } else {
    invMat <-
      opts$invTransform |>
      as.character() |>
      gsub("\\bx(\\d+)\\b", "x[,\\1]", x=_) |>
      gsub("\\bz(\\d+)\\b", "z[,\\1]", x=_) |>
      paste0(collapse=", ") |>
      paste0("cbind(", ...=_, ")") |>
      str2expression()
    out$matrix$invTransform <- function(x, z) eval(invMat)
    invVec <-
      opts$invTransform |>
      as.character() |>
      gsub("\\bx(\\d+)\\b", "x[\\1]", x=_) |>
      gsub("\\bz(\\d+)\\b", "z[\\1]", x=_) |>
      paste0(collapse=", ") |>
      paste0("c(", ...=_, ")") |>
      str2expression()
    out$vector$invTransform <- function(x, z) eval(invVec)
  }

  return(out)
}


prepareParmsGlobalLm <- function(parms, opts) {
  opts <- asOpts(opts, c("GlobalLm", "DerivFun"))
  lmFuns <- buildLmFuns(opts)
  z <- lmFuns$matrix$transform(parms$trajs$state, parms$trajs$deriv)
  coef <- lapply(seq_len(ncol(z)), \(j) {
    X <- lmFuns$matrix$features(parms$trajs$state, j)
    DEEButil::saveSolve(crossprod(X), crossprod(X, z[,j]))
  })
  parms$lmFuns <- lmFuns
  parms$coef <- coef
  return(parms)
}
