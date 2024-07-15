#' @export
print.colonization <- function(x, ...)
{
    cat("Minimum number of inferred colonization events",
        sep = "\n")
    cat(x$Total,
        sep = "\n")
}
