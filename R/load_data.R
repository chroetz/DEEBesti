#' @export
load_data <- function(file) {
  x <- readr::read_csv(file, col_types = readr::cols(), name_repair = "minimal")[-1]
  entries <- stringr::str_subset(names(x), "^u[0-9]+")
  x$u <- do.call(cbind, x[entries])
  x[entries] <- NULL
  entries <- stringr::str_subset(names(x), "^du[0-9]+")
  x$du <- do.call(cbind, x[entries])
  x[entries] <- NULL
  return(x)
}
