#' @export
summary.maxCol <- function(object) {
    cat("Final curve of colonization events in genetic estimator",
        sep = "\n")
    if(object$FormulaGen == "NA") {
        cat("Not possible to estimate the maximum of colonization events",
            sep = "\n")
    } else {
        cat(paste("Formula:",
                  object$FormulaGen,
                  sep = " "),
            sep = "\n")
        print(object$ParametersGen)
    }

    cat("",
        sep = "\n")
    cat("Final curve of colonization events in field estimator",
        sep = "\n")
    if (object$FormulaField == "NA") {
        cat("Not possible to estimate the maximum of colonization events",
            sep = "\n")
    } else {
        cat(paste("Formula:",
                  object$FormulaField,
                  sep = " "),
            sep = "\n")
        print(object$ParametersField)
    }
    cat("",
        sep = "\n")
    cat(paste("Interval of conficende of ",
              object$ConfintLevel * 100,
              "%",
              sep = ""),
        sep = "\n")
    cat("",
        sep = "\n")
    cat(paste("Deleted the",
              object$DeletedData * 100,
              "% of extreme values"),
        sep = "\n")
}
