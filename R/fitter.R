getFitter <- function(name) {
  switch(
    name,
    LocalConst = fitLocalConst,
    LocalLinear = fitLocalLinear,
    stop("Unknown fitter name ", name)
  )
}
