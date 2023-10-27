# Startup
#' .onAttach start message
#' @param libname defunct
#' @param pkgname defunct
#' @return invisible()
.onAttach <- function(libname, pkgname) {
	start_message <- "PAICE 1.0.0\n"
    packageStartupMessage(start_message)
    invisible()
}
