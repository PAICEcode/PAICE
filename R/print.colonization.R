#' @export
print.colonization <- function(object)
{
    cat("Total of inferred colonization events",
        sep = "\n")
    cat(object$Total,
        sep = "\n")
}
