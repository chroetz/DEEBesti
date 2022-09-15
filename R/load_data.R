#' @export
load_data <- function(file) {
  x <- readr::read_csv(file, col_types = readr::cols(), name_repair = "minimal")[-1]
  entries <- stringr::str_subset(names(x), "^state[0-9]+")
  x$state <- do.call(cbind, x[entries])
  x[entries] <- NULL
  return(x)
}
