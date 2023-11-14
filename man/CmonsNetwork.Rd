\name{CmonsNetwork}

\alias{CmonsNetwork}

\docType{data}

\title{Genealogical relationship of \emph{Cistus monspeliensis} haplotypes}

\description{
    Genealogy of Canarian haplotypes of \emph{Cistus monspeliensis}.
}

\usage{
data("CmonsNetwork")
}

\format{
    A data frame containing the genalogy of \emph{Cistus monspeliensis}
    in the Canary Islands. Each row indicates the connection between each
    haplotype and its ancestral haplotype. The first column is the name of the
    haplotype, the second column is the name of its ancestral haplotype and the
    third column indicates the number of variable positions that change between
    both haplotypes. The ancestral haplotype in the archipelago (haplotype C1)
    is connected to the outgroup (\code{"OUT"}), and is be located in the first
    row of the genealogy.
}

\details{
    This dataset was taken from Coello et al. (2021). It was constructed using
    three ptDNA regions and 37 populations from the Canarian archipelago.
}

\references{
    Coello, A.J., Fernandez-Mazuecos, M., Garcia-Verdugo, C., Vargas, P.
    (2021).
    Phylogeographic sampling guided by species distribution modeling reveals the
    Quaternary history of the Mediterranean-Canarian \emph{Cistus monspeliensis}
    (Cistaceae).
    \emph{Journal of Systematics and Evolution, 59(2)}, 262-277.
    DOI: \href{https://onlinelibrary.wiley.com/doi/10.1111/jse.12570}{10.1111/jse.12570}
}
\examples{
data(CmonsNetwork)
CmonsNetwork # Show data frame
}

\keyword{ datasets }
