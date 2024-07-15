#' @export
print.maxCol <- function(x, ...) {
    cat("Asymptote estimators of colonization events:",
        sep = "\n")
    print(x$Summary)
    cat("",
        sep = "\n")
    cat(paste("Minimum and maximum determined by a confidence interval of ",
              x$ConfintLevel * 100,
              "%",
              sep = ""),
        sep = "\n")
}
