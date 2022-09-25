prepareMemory <- function(opts, n) {
  opts <- asOpts(opts, "Method")
  method <- getClassAt(opts, 2)
  if (method == "Altopi") utils::clrhash(altopiMemory)
}
