\name{geneticResampling}

\alias{geneticResampling}

\title{Simulate genetic sampling effort reduction}

\description{
    This function reduces the resolution of the genealogy by supressing a
    variable position in the genealogy. It simulates a lower level of genetic
    sampling.
}

\usage{
geneticResample(data, network, position)
}

\arguments{
    \item{data}{
        a data frame containing the occurrence matrix of haplotypes in the
        islands of an archipelago (applicable to any island-like system). The
        two first columns indicate islands and populations sampled. Successive
        columns indicate haplotype occurrences (one column per haplotype). If
        present, missing haplotypes must also be included (i.e. columns without
        occurrences).
    }
    \item{network}{
        a data frame containing the genealogy of haplotypes. The first column
        indicates the haplotype, the second column indicates its ancestral
        haplotype and the third column indicates the variable position changed
        between an haplotype and its ancestral haplotype. If present, missing haplotypes must also be
        included. The ancestral haplotype must be connected to an outgroup named
        \code{"OUT"}, located in the first row of the data frame, and has a
        variable position not shared with other connections.
    }
    \item{position}{
        numeric. Indicates the variable position that will be deleted in the
        simplified data.
    }
}

\details{
    To simulate a lower level of genetic sampling, this function deletes a
    variable position from the original data and thus simplifies the genealogy.
    \code{geneticResampling} generates a new occurrence matrix of haplotypes and
    a new genealogy without the variable position previously indicated and
    merging ancestral and derived haplotypes separated by this variable
    position. If more than one connection are defined by the variable position
    indicated, this function deletes all connections with this variable
    position. This function works for both observed and missing haplotypes.
}

\value{
    \code{geneticResampling} returns a list containing the new occurrence matrix
    of haplotypes and the new genealogy after deleting the variable position
    indicated. The returned object contains the following components:

    \item{data}{
        a data frame containing the new occurrence matrix of haplotypes after
        removing the variable position indicated.
    }

    \item{network}{
        a data frame containing the new genealogy after removing the variable
        position indicated.
    }
}

\note{
    If the variable position corresponds to the connection between the ancestral
    haplotype in the archipelago and the outgroup (denoted as \code{"OUT"}), no
    change is effected as the ancestral haplotype stays connected to the
    outgroup.

    This function works inside \code{\link{rarecol}}.
}

\seealso{
    \code{\link{rarecol}} to build a rarefaction curve of colonization events.
}

\examples{
data(CmonsData)
data(CmonsNetwork)
# Delete position 462 of Cistus monspeliensis data
newdata <- geneticResampling(CmonsData, CmonsNetwork, 462)
newdata$data # New presences matrix of haplotypes
newdata$network # New genealogy
}

\keyword{ datagen }
