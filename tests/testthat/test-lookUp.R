findIdxDirect <- function(query, state, h) {
  dstAll <- distToVec(state, query)
  selAll <- dstAll <= h
  sort(which(selAll))
}

testRandomLookUp <- function(n, d, h, reps) {
  trajs <- makeDerivTrajs(matrix(rnorm(n*d), ncol=d))
  look <- buildLookUpGrid(trajs, distance=h)
  for (i in seq_len(reps)) {
    query <- rnorm(d)
    expect_identical(
      sort(lookAround(look, query, h)),
      findIdxDirect(query, trajs$state, h))
  }
}

test_that("random", {
  for (d in 1:4) for (n in 10^(0:3)) testRandomLookUp(n, d, 0.1, 20)
})
