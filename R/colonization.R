#' @export
#' @importFrom stats aggregate
colonization <- function(data, network) {
    # Colapse data by island and eliminate the population column
    data_by_island <- aggregate(x   = data[, -c(1, 2)],
                                by  = list(data[[1]]),
                                FUN = sum)
    names(data_by_island)[1] <- names(data)[1]

    # Transform colapsed data to presence/absence
    data_by_island[, -1] <- (data_by_island[, -1] >= 1) * 1

    # Calculate haplotypes occurred in each island
    if (is.data.frame(data_by_island[, -1]) == TRUE) {
        total_per_island <- apply(X      = data_by_island[, -1],
                                  MARGIN = 1,
                                  FUN    = sum)
    } else {
        # If data_by_island is not a data frame (only one haplotype exits) do not
        # colpase data
        total_per_island <- data_by_island[, -1]
    }
    
    # Calculate islands occurred
    total_per_island <- (total_per_island >= 1) * 1

    if (sum(total_per_island) != 1) {
        # Check the number of haplotypes in the data set
        total_per_hap <- apply(X      = as.data.frame(data_by_island[, -1]),
                               MARGIN = 2,
                               FUN    = sum)
        total_per_hap <- (total_per_hap >= 1) * 1

        if (sum(total_per_hap) != 1) {
            # If there are more than one haplotype in the data set

            # Continue preparing data

            # Delete missing haplotypes
            total <- colSums(data_by_island[, -1])
            hap_data <- cbind(data_by_island[, 1],
                              data_by_island[, names(total)[total >= 1]])
            names(hap_data)[1] <- names(data_by_island)[1]

            # Create a object with missing haplotypes name
            missing <- names(total)[total == 0]

            # Calculate whithin colonization events

            # Calculate total of islands in which each haplotype occurs
            island_hap <- colSums(hap_data[, -1])

            # Calculate whithin colonization events
            within <- island_hap - 1

            # Add an additional value in the begging for the outgroup
            within <- c(0, within)
            names(within)[1] <- "OUT"

            # Transform data for colonization events between

            # Search real ancestral haplotype to each haplotype

            # Create objects for loop
            haplotype <- c()
            ancestral <- c()

            for (h in 2:length(names(hap_data))) {
                # Determine haplotype
                haplotype2 <- names(hap_data)[h]

                # Prepare haplotype for loop
                haplotype3 <- haplotype2

                repeat {
                    # Determine ancestral haplotype name
                    ancestral2 <- as.vector(network[network[, 1] == haplotype3,
                                                    2])

                    # Check if the ancestral haplotype is the outgroup
                    # If it is the outgroup stop the loop
                    if (length(names(total)[names(total) == ancestral2]) == 0) {
                              break
                    } else {
                        # Determine number of islands occurred by the ancestral
                        # haplotype (if it is 0 then it is a missing haplotype)
                        ancestral_islands <- total[names(total) == ancestral2]
                        ancestral_islands <- as.vector(ancestral_islands)

                        if (ancestral_islands >= 1) {
                            break
                        } else {
                            # In case of missing haplotype as ancestral repeat
                            # loop with missing haplotype as haplotype in
                            # question
                            haplotype3 <- ancestral2
                        }
                    }
                }

                # Save the result
                haplotype <- c(haplotype, haplotype2)
                ancestral <- c(ancestral, ancestral2)
            }

            # Create a summary data frame for ancestry
            real_ancestry <- data.frame(haplotype, ancestral)

            # Prepare an object for loop
            sharing <- c()

            # Check if the haplotype share island with its ancestral haplotype
            for (h in 1:length(real_ancestry$haplotype)) {
                # Search haplotype
                haplotype <- real_ancestry$haplotype[h]

                # Search ancestrla haplotype
                ancestral <- real_ancestry$ancestral[h]

                # If ancestral is the outgroup do not compare
                if(ancestral != "OUT") {
                    # Search occurrences
                    haplotype_occ <- hap_data[, names(hap_data) == haplotype]
                    ancestral_occ <- hap_data[, names(hap_data) == ancestral]

                    # Compare islands
                    comparison <- sum(haplotype_occ * ancestral_occ)

                    # Indicate if haplotype share islands with its ancestral
                        if(comparison == 0) {
                            sharing <- c(sharing, "No")
                        } else {
                            sharing <- c(sharing, "Yes")
                        }
                } else {
                    sharing <- c(sharing, "OUT")
                }
            }

            # Add this information to the data frame
            real_ancestry <- cbind(real_ancestry, sharing)

            # Calculate derived colonization events
            # Prepare objects for loop
            derived <- c()
            marked <- c()

            for (h in 2:length(names(hap_data))) {
                haplotype <- names(hap_data)[h]

                # Check if this haplotype shares islands with its ancestral
                sharing <- real_ancestry$sharing[real_ancestry$haplotype ==
                                                 haplotype]
                sharing <- as.character(sharing)

                if (sharing == "Yes") {
                    # If the haplotype shares islands with its ancestral
                    # derived colonization events are not possible
                    derived <- c(derived, 0)
                } else if (sharing == "OUT") {
                    # If the ancestral haplotype is the outgroup check if the
                    # outgroup has more than one derived haplotypes
                    d <- real_ancestry$haplotype[real_ancestry$ancestral ==
                                                 "OUT"]
                    derived_out <- as.character(d)

                    # There are not derived colonizations from outgroup
                    # (only it is possible to infer ancestral from the outgroup)
                    derived <- c(derived, 0)

                    if (length(derived_out) >= 2) {
                        # Mark outgroup if it is not marked yet to calculate
                        # below ancestral colonizations
                        marked <- c(marked, "OUT")
                    }
                } else {
                    # Check the ancestral haplotype
                    anc <- real_ancestry$ancestral[real_ancestry$haplotype ==
                                                   haplotype]
                    ancestral <- as.character(anc)

                    # Check derived haplotypes from the ancestral haplotype
                    derh <- real_ancestry$haplotype[real_ancestry$ancestral ==
                                                    ancestral]
                    derived_haps <- as.character(derh)

                    if (length(derived_haps) == 1) {
# If the ancestral haplotype only shares islands with
# the haplotype infer one derived colonization event
                              derived <- c(derived, 1)
                    } else {
                        # Prepare an object for loop
                        delete <- c()

                        # If the ancestral haplotype derives in more than one
                        # haplotype, delete derived haplotypes that shares
                        # islands with the ancestral
                        for (i in 1:length(derived_haps)) {
                            # Search derived haplotype
                            ra <- real_ancestry
                            derived_sharing <- ra$sharing[ra$haplotype ==
                                                          derived_haps[i]]
                            derived_sharing <- as.character(derived_sharing)

                            # If derived haplotype shares with its ancestral
                            # delete it
                            if (derived_sharing == "Yes") {
                                delete <- c(delete,
                                            as.character(derived_haps[i]))
                            }
                        }

                        # Delete haplotypes marked
                        for (i in 1:length(delete)) {
                            dhcar <- as.character(derived_haps)
                            derived_haps <- derived_haps[dhcar != delete[i]]
                        }

                        if (length(derived_haps) == 1) {
                            # If after deleted shared derived haplotypes there
                            # are only one derived haplotype, infer one derived
                            # colonization event
                            derived <- c(derived, 1)
                        } else {
                            # Mark the ancestral haplotype if it is not market
                            # yet to calculate below ancestral colonization
                            # events and infer 0 derived colonization events
                            marked <- c(marked, ancestral)
                            derived <- c(derived, 0)
                        }
                    }
                }
            }

            # Add an additional value for the outgroup
            derived <- c(0, derived)

            # Create a data frame to unify colonization events already
            # calculated
            infer_col <- cbind(within, derived)

            # Infer ancestral colonization events
            
            # Calculate ancestral colonizations if there are marked ancestral
            # colonizators

            # Prepare an object for loop
            ancestral <- rep(x     = 0,
                             times = length(row.names(infer_col)))
            names(ancestral) <- row.names(infer_col)

            if (length(marked) >= 1) {
                # Determine ancestral haplotypes that produces ancestral
                # colonization events
                ancestral_colonizator <- unique(marked)

                for (h in 1:length(ancestral_colonizator)) {
                    # Extract ancestral colonizator occurrences
                    ancestral_occ <- hap_data[, names(hap_data) ==
                                                ancestral_colonizator[h]]

                    # When the ancestral colonizator is the outgroup subtract
                    # one ancestral colonization to do not consider arrival to
                    # the archipelago
                    if (ancestral_colonizator[h] == "OUT") {
                        ancestral[1] <- -1
                    }

                    # Search derived haplotypes
                    ra <- real_ancestry
                    derived_haps <- ra$haplotype[ra$ancestral ==
                                                 ancestral_colonizator[h]]
                    derived_haps <- as.character(derived_haps)

                    # Delete haplotypes that share islands with the ancestral
                    # haplotype (as above)
                    
                    # Prepare an object for loop
                    delete <- c()

                    # If the ancestral haplotype derives in more than one
                    # haplotype, delete derived haplotypes that shares islands
                    # with the ancestral
                    for (i in 1:length(derived_haps)) {
                        # Search derived haplotype
                        ra <- real_ancestry
                        derived_sharing <- ra$sharing[ra$haplotype ==
                                                      derived_haps[i]]
                        derived_sharing <- as.character(derived_sharing)

                        # If derived haplotype shares with its ancestral delete
                        # it
                        if (derived_sharing == "Yes") {
                            delete <- c(delete, as.character(derived_haps[i]))
                        }
                    }
                    
                    # Only delete haplotypes if there are haplotypes to delete
                    if (length(delete) >= 1) {
                        # Delete haplotypes marked
                        for (i in 1:length(delete)) {
                            dh <- as.character(derived_haps)
                            derived_haps <- derived_haps[dh != delete[i]]
                        }
                    }

                    # Search occurrences for derived haplotypes
                    for (i in 1:length(derived_haps)) {
                        x <- hap_data[, names(hap_data) == derived_haps[i]]

                        if (i == 1) {
                            derived_occ <- x
                        } else {
                            derived_occ <- cbind(derived_occ, x)
                        }
                    }
                    derived_occ <- as.data.frame(derived_occ)
                    row.names(derived_occ) <- as.vector(data_by_island[[1]])
                    names(derived_occ) <- derived_haps

                    # Object to stop the loop
                    stop <- "No"

                    repeat {
                        if (stop == "Yes") {
                            break
                        }

                        # Determine the number of derived haplotype per island
                        if (is.data.frame(derived_occ) == FALSE) {
                            # When there are only one derived haplotype, add the
                            # las ancestral colonization and stop the loop
                            ancestral[names(ancestral) ==
                                      ancestral_colonizator[h]] <-
                                ancestral[names(ancestral) ==
                                          ancestral_colonizator[h]] + 1
                            stop <- "Yes"
                        } else {
                            islands <- apply(X      = derived_occ,
                                             MARGIN = 1,
                                             FUN    = sum)

                            # Extract island names
                            names(islands) <- as.vector(data_by_island[[1]])

                            # Determine the maximum number of haplotypes
                            # occurred together in one island
                            max_occ <- max(islands)

                            if (max_occ == 0) {
                                # Stop when there are not more islands with
                                # haplotypes
                                stop <- "Yes"
                            } else {
                                # Determine islands with maximum number of
                                # different haplotypes
                                isl_max <- names(islands)[islands == max_occ]

                                # Select one random island that host the maximum
                                # number of different haplotypes (when there are
                                # only one islands this step has not effect)
                                isl_max <- sample(x    = isl_max,
                                                  size = 1)

                                # Determine haplotypes present in this island
                                # and delete it
                                doc <- row.names(derived_occ)
                                deleted_haps <- derived_occ[doc == isl_max, ]

                                # Prepare an object for loop
                                deleted_haps_erase <- c()

                                # Mark haplotypes to deleted
                                for (j in 1:length(names(deleted_haps))) {
                                    hap <- deleted_haps[, j]

                                    if (hap >= 1) {
                                        deleted_haps_erase <-
                                            c(deleted_haps_erase,
                                              names(deleted_haps[j]))
                                    }
                                }

                                # Delete haplotypes marked
                                for (j in 1:length(deleted_haps_erase)) {
                                    if (is.data.frame(derived_occ) == TRUE) {
                                        derived_occ <-
                                            derived_occ[, names(derived_occ) !=
                                                        deleted_haps_erase[j]]
                                    } else {
                                        # If there are only one haplotype
                                        # transform derived_occ, add the last
                                        # ancestral colonization and stop loop
                                        prl <-
                                            ancestral[names(ancestral) ==
                                                      ancestral_colonizator[h]]
                                        ancestral[names(ancestral) ==
                                                  ancestral_colonizator[h]] <-
                                            prl + 1
                                        stop <- "Yes"
                                    }
                                }

                                if (stop != "Yes"){
                                    # Add one ancestral colonization
                                    ancestral[names(ancestral) ==
                                              ancestral_colonizator[h]] <-
                                        ancestral[names(ancestral) ==
                                                  ancestral_colonizator[h]] + 1
                                }
                            }
                        }
                    }
                }
            }

            # Add information to the data frame of inferred colonization events
            infer_col <- cbind(infer_col, ancestral)
        } else {
            # If there are only one haplotype in the data set it can only be
            # possible to infer whithin colonization events

            # Continue preparing data

            # Delete missing haplotypes
            total <- colSums(as.data.frame(data_by_island[, -1]))
            names(total) <- names(data_by_island)[-1]

            hap_data <- as.data.frame(cbind(as.character(data_by_island[, 1]),
                                      data_by_island[, names(total)[total >=
                                                     1]]))
            names(hap_data)[1] <- names(data_by_island)[1]
            names(hap_data)[2] <- names(total)[total >= 1]

            # Calculate whithin colonization events

            # Calculate total of islands which each haplotype occurred
            island_hap <- sum(as.numeric(as.vector(hap_data[, -1])))

            # Calculate within colonization events
            within <- island_hap - 1
            names(within) <- names(total)[total >= 1]

            # Add an additional value in the begging for the outgroup
            within <- c(0, within)
            names(within)[1] <- "OUT"

            # There are not colonizations between (derived and ancestral)
            derived <- c(0, 0)
            ancestral <- c(0, 0)
            infer_col <- data.frame(within, derived, ancestral)
        }
    } else {
        within <- 0
        derived <- 0
        ancestral <- 0
        infer_col <- data.frame(within, derived, ancestral)
    }
    
    # Summary data

    # Calculate the number of islands in the archipelago
    islands <- data[apply(X      = as.data.frame(data[, -c(1, 2)]),
                          MARGIN = 1,
                          FUN    = sum) >= 1,
                    1]
    islands <- length(as.vector(unique(islands)))

     # Calcualte the number of populations used
     populations <- length(data[apply(X      = as.data.frame(data[, -c(1, 2)]),
                                      MARGIN = 1,
                                      FUN    = sum) >= 1,
                                2])

     # Calculate the number of individuals used
     individuals <- sum(data[apply(X      = as.data.frame(data[, -c(1, 2)]),
                                   MARGIN = 1,
                                   FUN    = sum) >= 1,
                             -c(1, 2)])

     # Calculate the number of haplotypes used
     hap_ind <- apply(X      = as.data.frame(data[, -c(1, 2)]),
                      MARGIN = 2,
                      FUN    = sum)
     hap_ind <- hap_ind[hap_ind >= 1]
     haplotypes <- length(hap_ind)

     summarydata <- matrix(data     = c(islands, populations, individuals,
                                        haplotypes),
                           nrow     = 4,
                           ncol     = 1,
                           byrow    = FALSE,
                           dimnames = list(c("Islands", "Populations",
                                             "Individuals", "Haplotypes"),
                                           c("Total")))

     # Format data about colonizations
     infer_col <- as.data.frame(infer_col)
     names(infer_col) <- c("c1", "c2", "c3")
     colonizationsbyhaplotype <- infer_col

     colonizationsbycomponent <- apply(X      = infer_col,
                                       MARGIN = 2,
                                       FUN    = sum)

     totalcolonizations <- matrix(data     = sum(colonizationsbycomponent),
                                  nrow     = 1,
                                  ncol     = 1,
                                  byrow    = FALSE,
                                  dimnames = list(c("Total of colonizations:"),
                                                  c("")))

     # Group result in a single object
     result <- list(Summary               = summarydata,
                    ColonizationHaplotype = colonizationsbyhaplotype,
                    ColonizationComponent = colonizationsbycomponent,
                    Total                 = totalcolonizations)
     class(result) <- "colonization"
     return(result)
}
