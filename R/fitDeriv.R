buildFitter <- function(opts) {
  opts <- asOpts(opts, "FitDeriv")
  name <- getClassAt(opts, 2)
  switch(
    name,
    LocalConst = \(x, y) fitterLocalConst(x, y, opts$bandwidth, getKernel(opts$kernel)),
    LocalLinear = \(x, y) fitterLocalLinear(x, y, opts$bandwidth, getKernel(opts$kernel)),
    GaussianProcess = \(x, y) fitterGaussianProcess(
      x, y, opts$bandwidth, getKernel(opts$kernel), opts$regulation, opts$neighbors),
    stop("Unknown fitter name ", name)
  )
}
