#' @export
#' @importFrom stats aggregate
#' @importFrom stats confint
#' @importFrom stats nls
#' @importFrom stats quantile
maxCol <- function(data, level = 0.95, del = 0.05, method = 1) {

    # Check parameters
    if (level >= 1 | level < 0) {
        stop("'level' must be a number between 0 and 1")
    }
    if (del >= 1 | del < 0) {
        stop("'del' must be a number between 0 and 1")
    }
    if(method != 1) {
        if(method != 0) {
            stop("'method' must be 1 or 0")
        }
    }

    # Create a function to delete extreme values
    fun <- function(x, del) {
        # Minimum and maximum values of the vector
        n <- quantile(x,
                      probs = c(del / 2,
                                1 - del / 2))
        # Filter data
        x <- x[x >= n[1] & x <= n[2]]

        # Calculate mean value
        result <- mean(x)
        return(result)
    }

    # First fit curve in field and then in genome

    # Assign data of genetic replicates as final curve
    data_gen <- data[[1]]

    g_replicate <- unique(data_gen$GeneticReplicate)

    # Determine number of variable positions
    n_pos <- unique(data_gen$VariablePositions)
    n_pos <- sort(n_pos)

    # Create objects to save data
    col_max <- vector(mode = "numeric")
    pos <- vector(mode = "numeric")

    for (g in g_replicate) {
        data_gen2 <- data_gen[data_gen$GeneticReplicate == g, ]
        for (n in n_pos) {
            data_gen3 <- data_gen2[data_gen2$VariablePositions == n, ]
            data_gen3 <- aggregate(x   = data_gen3,
                                   by  = list(data_gen3$Populations),
                                   FUN = mean)
            result <- tryCatch({
                    eq <- TotalColonizationEvents ~ col_max *
                          (Populations - 1) / (K + Populations - 1)
                    cmax <- max(data_gen3$TotalColonizationEvents)
                    kmax <- max(data_gen3$TotalColonizationEvents) / 2
                    model <- nls(formula = eq,
                                 data    = data_gen3,
                                 start   = list(col_max = cmax,
                                                K       = kmax))
                    result <- c(summary(model)$coefficients[1, 1],
                                summary(model)$coefficients[2, 1])
                }, error = function(err) {
                    result <- NA
                })

            col_max <- c(col_max, result[1])
            pos <- c(pos, n)
        }
    }

    MM_data_gen <- data.frame(col_max, pos)
    MM_data_gen <- MM_data_gen[is.na(MM_data_gen$col_max) == FALSE, ]
    
    # Delete the extreme values
    MM_data_gen <- aggregate(x   = MM_data_gen,
                             by  = list(MM_data_gen$pos),
                             FUN = function(x) fun(x, del))

    result <- tryCatch({
            MM_nls_gen <- nls(formula = col_max ~ pos * M / (K + pos) + c,
                              data    = MM_data_gen,
                              start   = list(M = max(MM_data_gen$col_max),
                                             K = max(MM_data_gen$col_max) / 2,
                                             c = 0))
            mean_gen <- summary(MM_nls_gen)$coefficients
            CI_gen <- confint(object = MM_nls_gen,
                              level  = level)
            M_gen <- mean_gen[1, 1]
            K_gen <- mean_gen[2, 1]
            c_gen <- mean_gen[3, 1]
            M_gen_min <- CI_gen[1, 1]
            M_gen_max <- CI_gen[1, 2]
            K_gen_min <- CI_gen[2, 1]
            K_gen_max <- CI_gen[2, 2]
            c_gen_min <- CI_gen[3, 1]
            c_gen_max <- CI_gen[3, 2]
            col_max_total <- M_gen + c_gen
            col_max_max <- M_gen_max + c_gen_max
            col_max_min <- M_gen_min + c_gen_min

            MM_nls_gen <- matrix(data     = c(col_max_total,
                                              col_max_min, col_max_max),
                                 nrow     = 1,
                                 ncol     = 3,
                                 byrow    = TRUE,
                                 dimnames = list(c("Genetic estimator"),
                                                 c("Mean", "Min", "Max")))

            Gen_params <- matrix(data     = c(M_gen, M_gen_min, M_gen_max,
                                              K_gen, K_gen_min, K_gen_max,
                                              c_gen, c_gen_min, c_gen_max),
                                 nrow     = 3,
                                 ncol     = 3,
                                 byrow    = FALSE,
                                 dimnames = list(c("mean", "min", "max"),
                                                 c("m", "k", "c0")))
            formula_gen <- "c = m * positons / (k + positions) + c0"
            result <- list(MM_nls_gen, Gen_params, formula_gen)
        }, error = function(err) {
            result <- tryCatch({
                if (method == 1) {
                    Mmax <- max(MM_data_gen$col_max)
                    Kmax <- max(MM_data_gen$col_max) / 2
                    MM_nls_gen <- nls(formula = col_max ~ pos * M / (K + pos) +
                                                MM_data_gen$col_max[1],
                                      data    = MM_data_gen,
                                      start   = list(M = Mmax,
                                                     K = Kmax))
                    mean_gen <- summary(MM_nls_gen)$coefficients
                    CI_gen <- confint(object = MM_nls_gen,
                                      level  = level)
                    M_gen <- mean_gen[1, 1]
                    K_gen <- mean_gen[2, 1]
                    c_gen <- MM_data_gen$col_max[1]
                    M_gen_min <- CI_gen[1, 1]
                    M_gen_max <- CI_gen[1, 2]
                    K_gen_min <- CI_gen[2, 1]
                    K_gen_max <- CI_gen[2, 2]
                    c_gen_min <- NA
                    c_gen_max <- NA
                    col_max_total <- M_gen + c_gen
                    col_max_max <- M_gen_max + c_gen
                    col_max_min <- M_gen_min + c_gen

                    MM_nls_gen <- matrix(data     = c(col_max_total,
                                                      col_max_min, col_max_max),
                                         nrow     = 1,
                                         ncol     = 3,
                                         byrow    = TRUE,
                                         dimnames = list(c("Genetic estimator"),
                                         c("Mean", "Min", "Max")))
                    Gen_params <- matrix(data     = c(M_gen, M_gen_min,
                                                      M_gen_max,
                                                      K_gen, K_gen_min,
                                                      K_gen_max,
                                                      c_gen, c_gen_min,
                                                      c_gen_max),
                                         nrow     = 3,
                                         ncol     = 3,
                                         byrow    = FALSE,
                                         dimnames = list(c("mean", "min",
                                                           "max"),
                                                         c("m", "k", "c0")))
                    formula_gen <- "c = m * positons / (k + positions) + c0"
                    result <- list(MM_nls_gen, Gen_params, formula_gen)
                } else {
                    stop("Not calculate this case")
                }
            }, error = function(err) {
                text <- "Not possible to fit the curve of colonization events"
                result <- list(matrix(data = NA,
                                      nrow = 1,
                                      ncol = 3,
                                      dimnames = list(c("Genetic estimator"),
                                                      c("Mean", "Min", "Max"))),
                               text,
                               "NA")
            })
        })

    MM_data_gen <- MM_data_gen[, -1]
    names(MM_data_gen) <- c("colonization", "positions")
    MM_nls_gen <- result[[1]]
    Gen_params <- result[[2]]
    formula_gen <- result[[3]]

    # Then fit curve in field and then in genome

    # Assign data of field replicates as final curve
    data_field <- data[[2]]

    f_replicate <- unique(data_field$FieldReplicate)

    # Determine number of variable positions
    n_pop <- unique(data_field$Populations)
    n_pop <- sort(n_pop)

    # Create objects to save data
    col_max <- vector(mode = "numeric")
    pop <- vector(mode = "numeric")
    type <- vector(mode = "character")

    for (f in f_replicate) {
        data_field2 <- data_field[data_field$FieldReplicate == f, ]

        for (n in n_pop) {
            data_field3 <- data_field2[data_field2$Populations == n, ]
            data_field3 <- aggregate(x   = data_field3,
                                     by  = list(data_field3$VariablePositions),
                                     FUN = mean)

            if (min(data_field3$TotalColonizationEvents) ==
                max(data_field3$TotalColonizationEvents)) {
                
                c <- max(data_field3$TotalColonizationEvents)
                t <- "Fixed value"
                result <- list(c, t)
            } else {
                result <- tryCatch({
                        colmax <- max(data_field3$TotalColonizationEvents)
                        Kmax <- max(data_field3$TotalColonizationEvents) / 2
                        model <- nls(formula = TotalColonizationEvents ~
                                               (col_max * VariablePositions) /
                                               (K + VariablePositions) + c,
                                     data    = data_field3,
                                     start   = list(col_max = colmax,
                                                    K       = Kmax,
                                                    c       = 0))
                        c <- summary(model)$coefficients[1, 1] +
                             summary(model)$coefficients[3, 1]
                        t <- "Full model"
                        result <- list(c, t)
                    }, error = function(err) {
                        result <- tryCatch({
                                eq <- TotalColonizationEvents ~
                                      (col_max * VariablePositions) /
                                      (K + VariablePositions) +
                                      data_field3$TotalColonizationEvents[1]
                                cm <- max(data_field3$TotalColonizationEvents)
                                Km <- max(data_field3$TotalColonizationEvents) /
                                      2
                                if (method == 1) {
                                    model <- nls(formula = eq,
                                                 data    = data_field3,
                                                 start   = list(col_max = cm,
                                                                K       = Km))
                                c <- summary(model)$coefficients[1, 1] +
                                     data_field3$TotalColonizationEvents[1]
                                t <- "Not estimated intercept"
                                result <- list(c, t)
                            } else {
                                stop("Not calculate this case")
                            }
                        }, error = function(err) {
                            c <- NA
                            t <- NA
                            result <- list(c, t)
                    })
                })
            }
            col_max <- c(col_max, result[[1]])
            type <- c(type, result[[2]])
            pop <- c(pop, n)
        }
    }

    # After that infer curve in the field
    MM_data_field <- data.frame(col_max, pop, type)
    MM_data_field <- MM_data_field[is.na(MM_data_field$col_max) == FALSE, ]
    MM_data_field <- aggregate(x   = MM_data_field[, c(1, 2)],
                               by  = list(MM_data_field$pop),
                               FUN = function(x) fun(x, del))

    result <- tryCatch({
            MM_nls_field <- nls(formula = col_max ~ (pop - 1) * M /
                                          (K + pop - 1),
                                data    = MM_data_field,
                                start   = list(M = max(MM_data_field$col_max),
                                               K = max(MM_data_field$col_max) /
                                                   2))

            CI <- confint(object = MM_nls_field,
                          level  = level)

            M_field <- summary(MM_nls_field)$coefficients[1, 1]
            M_field_min <- CI[1, 1]
            M_field_max <- CI[1, 2]
            K_field <- summary(MM_nls_field)$coefficients[2, 1]
            K_field_min <- CI[2, 1]
            K_field_max <- CI[2, 2]
            MM_nls_field <- matrix(data     = c(M_field, M_field_min,
                                              M_field_max),
                                   nrow     = 1,
                                   ncol     = 3,
                                   byrow    = TRUE,
                                   dimnames = list(c("Field estimator"),
                                                   c("Mean", "Min", "Max")))
            Field_params <- matrix(data     = c(M_field, M_field_min,
                                                M_field_max, K_field,
                                                K_field_min, K_field_max),
                                   nrow     = 3,
                                   ncol     = 2,
                                   byrow    = FALSE,
                                   dimnames = list(c("mean", "min", "max"),
                                                   c("m", "k")))
            if (max(data[[1]]$Individuals) == max(data[[1]]$Populations)) {
                text <- "c = m * (individuals - 1) / (k + individuals - 1)"
                formula_field <- text
            } else {
                text <- "c = m * (populations - 1) / (k + populations - 1)"
                formula_field <- text
            }
            result <- list(MM_nls_field, Field_params, formula_field)
        }, error = function(err) {
            result <- list(matrix(data     = NA,
                                  nrow     = 1,
                                  ncol     = 3,
                                  byrow    = FALSE,
                                  dimnames = list(c("Field estimator"),
                                                  c("Mean", "Min", "Max"))),
                        "Not possible to fit the curve of colonization events",
                        "NA")
        })

    Field_params <- result[[2]]
    MM_summary <- rbind(MM_nls_gen, result[[1]])
    MM_data_field <- MM_data_field[, -1]
    if (max(data[[1]]$Individuals) == max(data[[1]]$Populations)) {
        names(MM_data_field) <- c("colonization", "individuals")
    } else {
        names(MM_data_field) <- c("colonization", "populations")
    }
    formula_field <- result[[3]]

    # Output

    result <- list(DataGen         = MM_data_gen,
                   FormulaGen      = formula_gen,
                   DataField       = MM_data_field,
                   FormulaField    = formula_field,
                   Summary         = MM_summary,
                   ParametersGen   = t(Gen_params),
                   ParametersField = t(Field_params),
                   ConfintLevel    = level,
                   DeletedData     = del)

    class(result) <- "maxCol"

    # Return the result
    return(result)
}
