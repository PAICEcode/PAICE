\name{CmonsData}

\alias{CmonsData}

\docType{data}

\title{Occurrence matrix of \emph{Cistus monspeliensis} in the Canary Islands}

\description{
    Data of \emph{Cistus monspeliensis} prepared to be used as example for the
    \pkg{PAICE} package.
}

\usage{
data(CmonsData)
}

\format{
    A data frame containing a presence matrix of \emph{Cistus monspeliensis}
    haplotypes in the Canary Islands extracted from Coello et al. (2021). Each
    row indicates the number of individuals of each haplotype occurring in each
    population. The first column indicates the island, the second column
    indicates the population and successive columns correspond to haplotypes in
    the island system. Missing haplotypes are also included but without any
    presence (haplotypes \code{m1} and \code{m2}).
}

\details{
    Data containing occurrences of each haplotype of \emph{Cistus monspeliensis}
    found in the Canary Islands. Data were taken from Coello et al. (2021). This
    dataset was constructed using three ptDNA regions and 37 populations from
    the Canarian archipelago.
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
data(CmonsData)
CmonsData # Show data frame
}

\keyword{ datasets }
