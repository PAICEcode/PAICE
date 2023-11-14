#' @export
# Simplify haplotypes by deleting variable position
geneticResampling <- function(data, network, position) {
    # Determine conections to delete
    coincidences <- network[network[, 3] == position, ]

    # If no coincidences are present stop the function
    if (nrow(coincidences) == 0) {
        stop(paste("Position ",
                   position,
                   " is not a variable position. Genetic resample is not ",
                   "possible",
                   sep = ""))
    }

    # Create vectors to save the process
    delhap <- vector(mode = "character")
    delanc <- vector(mode = "character")

    # Repeat the loop of simply network for each coincidence (it is possible to
    # have more than one edge with the same position in a network)
    repeat {
        # Extract haplotype of first coincidence
        hap <- as.character(coincidences[1, 1])

        # Extract ancestral
        ancestral <- as.character(coincidences[1, 2])

        # Check if the ancestral haplotype is already deleted
        if (length(delhap) > 0) {
            repeat {
                if (sum(delhap == ancestral) >= 1) {
                    ancestral <- delanc[delhap == ancestral]
                } else {
                    break
                }
            }
        }

        if (ancestral == "OUT") {
            # If the ancestral is the outgroup then delete this coincidence
            coincidences <- coincidences[-1, ]

            if (nrow(coincidences) == 0) break
        } else {
            # Search descendants
            descendants <- as.character(network[network[[2]] == hap, 1])

            repeat {
                # If there are not descendants stop the loop
                if (length(descendants) == 0) break

                # Repeat the proccess with all descendants
                
                # Assign ancestral to first descendant
                network[network[[1]] == descendants[1], 2] <- ancestral

                # Delete first descendant
                descendants <- descendants[-1]

            }
        }

        # Delete haplotype in network
        network <- network[network[[1]] != hap, ]

        # Add data of haplotype to ancestral in 'data' matrix
        data[, names(data) == ancestral] <-
            data[, names(data) == ancestral] + data[, names(data) == hap]

        # Delete haplotype in 'data' matrix
        data <- data[, names(data) != hap]

        # Save the first coincidence to check it
        delhap <- c(delhap, as.character(coincidences[1, 1]))
        delanc <- c(delanc, as.character(coincidences[1, 2]))

        # Delete first coincidence
        coincidences <- coincidences[-1, ]

        # If no more coincidences are present stop the loop
        if(nrow(coincidences) == 0) break
    }

    # Save new data and network
    result <- list(data = data, network = network)
    return(result)
}
