#' @export
#' @importFrom utils read.csv
read.rarecol <- function(gen, field) {
    gen <- read.csv(file = gen)
    field <- read.csv(file = field)
    result <- list(finaldataGen = gen,
                   finaldataFie = field)
    class(result) <- "rarecol"
    return(result)
}
