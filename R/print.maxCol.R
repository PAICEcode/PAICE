#' @export
print.maxCol <- function(object) {
    cat("Maximum of estimated colonization events:",
        sep = "\n")
    print(object$Summary)
    cat("",
        sep = "\n")
    cat(paste("Minimum and maximum determined by a interval of conficende of ",
              object$ConfintLevel * 100,
              "%",
              sep = ""),
        sep = "\n")
}
