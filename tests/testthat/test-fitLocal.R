test_that("fitterLocalLinear on vectors", {
  expect_equal(
    fitterLocalLinear(1:2, 1:2, bw=2, kernel=kernParab),
    matrix(1:2, ncol=1))
})
