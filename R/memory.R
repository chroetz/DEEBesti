prepareMemory <- function(method, n) {
  if (method == "Altopi")
    altopiMemory <<- utils::hashtab("identical", n)
}
