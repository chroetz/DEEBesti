.onLoad <- function(libname, pkgname) {
  options(warnPartialMatchDollar = TRUE)
  ConfigOpts::addPackageToPathDefaults(
    system.file("defaultOpts", package=pkgname, lib.loc=libname))
  invisible(NULL)
}
