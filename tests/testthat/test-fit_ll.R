test_that("fit_ll on vectors", {
  expect_equal(
    fit_ll(1:2, 1:2, bw=2, kernel=kern_parab),
    matrix(1:2, ncol=1))
})
