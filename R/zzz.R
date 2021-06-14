# Startup
#' .onAttach start message
#' @param libname defunct
#' @param pkgname defunct
#' @return invisible()
.onAttach <- function(libname, pkgname) {
	start_message <- "PAICE loaded\n\n"
    packageStartupMessage(start_message)
    invisible()
}
