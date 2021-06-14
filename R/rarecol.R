#' @export
#' @importFrom utils write.table
rarecol<- function(data, network, replicates_field = 10,
                   replicates_genetic = 10, mode = c(1, 2), monitor = TRUE,
                   file = NULL) {
    # Check labels in individuals/populations
    if (length(unique(data[[2]])) < length(data[[2]])) {
        stop("Two or more individuals/populations have the same name")
    }

    # Duplicate original objects
    original_data <- data
    original_network <- network

    if (is.null(file) == FALSE) {
        fileGen <- paste(file,
                         "gen.csv",
                         sep = "_")
        fileField <- paste(file,
                           "field.csv",
                           sep = "_")
    } else {
        fileGen <- NULL
        fileField <- NULL
    }

    # Approach in genome and then in field
    if (sum(mode == 1)) {
        # Create an empty matrix to save each point of resampling to infer the
        # maximum of colonization events
        finaldataGen <- matrix(data = NA,
                               nrow = 0,
                               ncol = 5)
        colnames(finaldataGen) <- c("Populations",
                                    "Individuals",
                                    "TotalColonizationEvents",
                                    "VariablePositions",
                                    "GeneticReplicate")
    
        if (is.null(fileGen) == FALSE) {
            write.table(x       = finaldataGen,
                        file    = fileGen,
                        sep     = ",",
                        qmethod = "double")
        }
    
        # Accumulation curves data generation
        for(genetic in 1:replicates_genetic) {
            # Replicates for genetic sampling
            if(monitor == TRUE) {
                cat(paste("[",
                          format(Sys.time(),
                                 "%d/%m/%y, %X"),
                          "] ",
                          genetic,
                          " genetic replicates",
                          sep = ""),
                    sep = "\n")
            }
    
            # Go back to original data
            data <- original_data
            network <- original_network
    
            # Create an object to stop genetic resampling
            end <- "no"
    
            # Determine unique variable positions and randomize it
            variable_positions <- sample(unique(network[[3]]))
    
            # Start with the most complete genetic data and delete variable
            # positions until there is no genetic data
            repeat {
                for (field in 1:replicates_field) {
                    # Determine number of populations
                    n_pop <- length(data[[1]])
    
                    # Determine order to add populations
                    pop_order <- sample(x       = 1:n_pop,
                                        replace = FALSE)
    
                    # Prepare data for loop
                    Populations <- c()
                    Individuals <- c()
                    TotalColonizationEvents <- c()
                    VariablePositions <- c()
                    GeneticReplicate <- c()
    
                    for (pop in 1:n_pop) {
                        # Infer colonization events from initial population to
                        # population pop in pop_order
                        data_col <- data[pop_order[1:pop], ]
                        col_i <- colonization(data    = data_col,
                                               network = network)
    
                        # Save information
                        Populations <- c(Populations, col_i$Summary[2, 1])
                        Individuals <- c(Individuals, col_i$Summary[3, 1])
                        TotalColonizationEvents <- c(TotalColonizationEvents,
                                                     col_i$Total[1, 1])
                        VariablePositions <- c(VariablePositions,
                                               length(variable_positions))
                        GeneticReplicate <- c(GeneticReplicate, genetic)
                    }
    
                    # Create a data frame with the information of the
                    # accumulation curves of colonization events
                    Summary <- data.frame(Populations, Individuals,
                                          TotalColonizationEvents,
                                          VariablePositions, GeneticReplicate)
    
                    if (is.null(fileGen) == TRUE) {
                             # Add data to finaldata object
                             finaldataGen <- rbind(finaldataGen, Summary)
                    } else {
                             write.table(x         = Summary,
                                         file      = fileGen,
                                         sep       = ",",
                                         qmethod   = "double",
                                         append    = TRUE,
                                         col.names = FALSE,
                                         row.names = FALSE)
                    }
                }
                
                if(monitor == TRUE) {
                    cat(paste("     [",
                              format(Sys.time(),
                                     "%d/%m/%y, %X"),
                              "] ",
                              length(variable_positions),
                              " positions",
                              sep = ""),
                        sep = "\n")
                }
                
                # Delete first variable position in datasets (but not in the
                # last iteration with only chorology data)
                if (length(variable_positions) > 0){
                    pos_nd <- variable_positions[1]
                    new_data <- geneticResampling(data     = data,
                                                  network  = network,
                                                  position = pos_nd)
    
                    # Save new data
                    data <- new_data[[1]]
                    network <- new_data[[2]]
                }
    
                # Delete first position for future iterations
                variable_positions <- variable_positions[-1]
    
                # If there are not more variable positions add colonization
                # events by chorology and stop the loop
                if (length(variable_positions) == 0) {
                    if (end == "no") {
                        end <- "yes"
                    } else {
                        break
                    }
                }
            }
        }
    
        if (is.null(fileGen) == FALSE) {
            finaldataGen <- fileGen
            names(finaldataGen) <- "Data saved in "
        }
    }
    
    # Approach in field and then in genome
    if (sum(mode == 2)) {
        # Repeat the process with the inverse approach (in this case field and
        # then genetic)
    
        # Create an empty matrix to save each point of resampling to infer the
        # maximum of colonization events
        finaldataFie <- matrix(data = NA,
                               nrow = 0,
                               ncol = 5)
        colnames(finaldataFie) <- c("Populations",
                                    "Individuals",
                                    "TotalColonizationEvents",
                                    "VariablePositions",
                                    "FieldReplicate")
    
        if (is.null(fileField) == FALSE) {
            write.table(x       = finaldataFie,
                        file    = fileField,
                        sep     = ",",
                        qmethod = "double")
        }
    
        for (field in 1:replicates_field) {
            # Replicates for field sampling
            if(monitor == TRUE) {
                cat(paste("[",
                          format(Sys.time(),
                                 "%d/%m/%y, %X"),
                          "] ",
                          field,
                          " field replicates",
                          sep = ""),
                    sep = "\n")
            }
    
            # Go back to original data
            data <- original_data
            network <- original_network
    
            populations <- as.vector(sample(x = data[[2]]))
    
            # Start with the most complete case and subtract one population each
            # time
            repeat {
                for (genetic in 1:replicates_genetic) {
                    # Recover data of this field replicate
                    data2 <- data
                    network2 <- network
    
                    # Determine unique variable positions and randomize it
                    variable_positions <- sample(unique(network2[[3]]))
    
                    # Create an object to finish genetic resampling
                    end <- "no"
    
                    # Prepare data for loop
                    Populations <- c()
                    Individuals <- c()
                    TotalColonizationEvents <- c()
                    VariablePositions <- c()
                    FieldReplicate <- c()
    
                    # Start with maximum genetic data and subtract one position
                    # each time
                    repeat {
                        # Infer colonizations from initial population to
                        # population pop in pop_order
                        col_i <- colonization(data    = data2,
                                               network = network2)
    
                        # Save information
                        Populations <- c(Populations, col_i$Summary[2, 1])
                        Individuals <- c(Individuals, col_i$Summary[3, 1])
                        TotalColonizationEvents <- c(TotalColonizationEvents,
                                                     col_i$Total[1, 1])
                        VariablePositions <- c(VariablePositions,
                                               length(variable_positions))
                        FieldReplicate <- c(FieldReplicate, field)
    
                        if (length(variable_positions) > 0) {
                            resampos <- variable_positions[1]
                            new_data <- geneticResampling(data     = data2,
                                                          network  = network2,
                                                          position = resampos)
    
                            # Save new data
                            data2 <- new_data[[1]]
                            network2 <- new_data[[2]]
                        }
    
                        # Delete first position for future iterations
                        variable_positions <- variable_positions[-1]
    
                        # If there are not more variable positions add
                        # colonization events by chorology and stop the loop
                        if (length(variable_positions) == 0) {
                            if (end == "no") {
                                end <- "yes"
                            } else {
                                break
                            }
                        }
                    }
    
                    # Create a data frame with the information of the
                    # accumulation curves of colonization events
                    Summary <- data.frame(Populations,
                                          Individuals,
                                          TotalColonizationEvents,
                                          VariablePositions,
                                          FieldReplicate)
    
                    if (is.null(fileField) == TRUE) {
                        # Add data to finaldata object
                        finaldataFie <- rbind(finaldataFie, Summary)
                    } else {
                        write.table(x         = Summary,
                                    file      = fileField,
                                    sep       = ",",
                                    qmethod   = "double",
                                    append    = TRUE,
                                    col.names = FALSE,
                                    row.names = FALSE)
                    }
                }
    
                if (monitor == TRUE) {
                    cat(paste("     [",
                              format(Sys.time(),
                                     "%d/%m/%y, %X"),
                              "] ",
                              length(populations),
                              " individuals/populations",
                              sep = ""),
                        sep = "\n")
                }
    
                # Delete one population
                data <- data[data[[2]] != populations[1], ]
    
                # Delete first population in the random vector
                populations <- populations[-1]
    
                if (length(populations) == 0) break
            }
        }
        
        if (is.null(fileField) == FALSE) {
            finaldataFie <- fileField
            names(finaldataFie) <- "Data saved in "
        }
    }
    
    # Output
    
    if (prod(mode == c(1, 2))) {
        finaldata <- list(finaldataGen   = finaldataGen,
                          finaldataField = finaldataFie)
    } else if(mode == 1) {
            finaldata <- list(finaldataGen = finaldataGen)
    } else {
            finaldata <- list(finaldataField = finaldataFie)
    }
    
    class(finaldata) <- "rarecol"
    if (is.null(file)) {
             return(finaldata)
    }
}
