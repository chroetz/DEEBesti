buildFitter <- function(opts) {
  opts <- asOpts(opts, "Fitter")
  name <- getClassAt(opts, 2)
  switch(
    name,
    LocalConst = \(x, y) fitLocalConst(x, y, opts$bandwidth, getKernel(opts$kernel)),
    LocalLinear = \(x, y) fitLocalLinear(x, y, opts$bandwidth, getKernel(opts$kernel)),
    stop("Unknown fitter name ", name)
  )
}
