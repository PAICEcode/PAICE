#' @export
print.maxCol <- function(x, ...) {
    cat("Maximum of estimated colonization events:",
        sep = "\n")
    print(x$Summary)
    cat("",
        sep = "\n")
    cat(paste("Minimum and maximum determined by a interval of conficende of ",
              x$ConfintLevel * 100,
              "%",
              sep = ""),
        sep = "\n")
}
