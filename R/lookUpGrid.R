buildLookUpGrid <- function(trajs, distance) {
  d <- ncol(trajs$state)
  n <- nrow(trajs$state)

  ranges <- sapply(seq_len(d), \(j) {
    rng <- range(trajs$state[,j])
    if (diff(rng) < sqrt(.Machine$double.eps)) {
      rng <- mean(rng) + c(-1,1)*sqrt(.Machine$double.eps)/2
    }
    return(rng)
  })
  centers <- apply(ranges, 2, \(r) (r[2]+r[1]) / 2)

  s <- apply(ranges, 2, \(r) ceiling((r[2]-r[1]) / distance / 2) + 1)

  getKey <- function(x) {
    res <- round((x-centers)/distance) + s + 1
    res[res<=0] <- 1
    sel <- res > 2*s+1
    res[sel] <- 2*s[sel]+1
    as.integer(res)
  }
  roundedState <-
    apply(trajs$state, 1, getKey, simplify=FALSE) |>
    unlist() |>
    matrix(nrow = n, byrow=TRUE)
  table <- hashtab()
  for (i in seq_len(n)) {
    key <- roundedState[i, ]
    value <- c(gethash(table, key), i)
    sethash(table, key, value)
  }
  return(list(table = table, getKey = getKey, trajs = trajs, distance = distance))
}

lookUp <- function(look, query) {
  key <- look$getKey(query)
  surrounding <- getSurronding(key)

  # Try to be faster than
  # collection <- apply(surrounding, 2, \(key) gethash(look$table, key), simplify=FALSE)
  len <- ncol(surrounding)
  d <- nrow(surrounding)
  collection <- vector(mode = "list", length = len)
  idx <- -(d-1):0
  for (i in seq_len(len)) {
    collection[[i]] <- gethash(look$table, surrounding[i*d+idx])
  }
  return(unlist(collection))
}

getSurronding <- function(key) {
  d <- length(key)
  out <- matrix(0L, ncol = 3^d, nrow = d)
  z <- seq_len(3^d)-1
  for (k in seq_len(d)) {
    out[k, ] <- key[k] + as.integer((z %/% 3^(k-1)) %% 3 - 1)
  }
  out
}

lookAround <- function(look, query, kMax) {
  idx <- lookUp(look, query)
  candidates <- look$trajs$state[idx, , drop=FALSE]
  dst <- distToVec(candidates, query)
  sel <- dst <= look$distance
  if (sum(sel) > kMax) {
    dstSel <- dst[sel]
    selSel <- .Internal(rank(dstSel, length(dstSel), "min")) <= kMax
    sel[sel] <- selSel
  }
  return(as.integer(idx[sel]))
}
