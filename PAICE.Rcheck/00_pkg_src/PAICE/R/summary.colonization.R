#' @export
summary.colonization <- function(object, ...) {
    cat("Summary of data used:",
        sep = "\n")
    print(object$Summary)
    cat("",
        sep = "\n")
    cat("Inference of colonization events: c = c1 + c2 + c3",
        sep = "\n")
    cat("",
        sep = "\n")
    cat("Colonization events by haplotype:",
        sep = "\n")
    print(object$ColonizationHaplotype)
    cat("",
        sep = "\n")
    cat("Colonization events by type:",
        sep = "\n")
    print(object$ColonizationComponent)
    cat("",
        sep = "\n")
    cat("Total of colonization events inferred:",
        sep = "\n")
    cat(paste("c =",
              object$Total,
              "colonization events",
              sep = " "),
        sep = "\n")
}
