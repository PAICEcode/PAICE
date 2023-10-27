#' @export
print.colonization <- function(x, ...)
{
    cat("Total of inferred colonization events",
        sep = "\n")
    cat(x$Total,
        sep = "\n")
}
