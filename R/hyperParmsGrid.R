getOptsGridList <- function(optsClass, ...) { # TODO: move to ConfigOpts
  default <- getDefaultOpts(optsClass)
  args <- list(...)
  isSingleOptsInDefault <- sapply(default, \(x) isOpts(x) && !inheritsOptsClass(x, "List"))
  optsNamesInDefault <- names(default)[isSingleOptsInDefault]
  isArgsOptsList <- sapply(args, \(x) inheritsOptsClass(x, "List"))
  argsNamesOptsList <- names(args)[isArgsOptsList]
  expandOptsListNames <- intersect(argsNamesOptsList, optsNamesInDefault)
  for (nm in expandOptsListNames) {
    lst <- args[[nm]]$list
    lstNames <- names(lst)
    attributes(lst) <- NULL
    names(lst) <- lstNames
    args[[nm]] <- lst
  }
  grid <- tidyr::expand_grid(!!!args)
  lst <- dfAsListOfOpts(grid, optsClass)
  makeOpts("List", list = lst)
}

dfAsListOfOpts <- function(df, optsClass) {
  lapply(seq_len(nrow(df)), \(i) {
    row <- getRow(df, i)
    opts <- do.call(makeOpts, c(list(optsClass), row))
    return(opts)
  })
}

getRow <- function(df, i) {
  row <- lapply(seq_len(ncol(df)), \(j) {
    v <- df[[i,j]]
    if (is.list(v)) {
      stopifnot(length(v) == 1)
      return(v[[1]])
    }
    v
  })
  names(row) <- names(df)
  return(row)
}
