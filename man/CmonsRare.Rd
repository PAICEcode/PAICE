\name{CmonsRare}

\alias{CmonsRare}

\docType{data}

\title{Simulated rarefaction curves of \emph{Cistus monspeliensis}}

\description{
    Simulated rarefaction curves to be used as example data for estimation of
    colonization events.
}

\usage{
data(CmonsRare)
}

\format{
    A list containing data of both genetic and field rarefaction curves. The
    first element corresponds to the genetic estimation and the second element
    corresponds to the field estimation.
}

\details{
    This dataset was constructed from \code{\link{CmonsData}} and
    \code{\link{CmonsNetwork}} with the following code:

    \preformatted{
set.seed(31)

CmonsRare <- rarecol(data = CmonsData, network = CmonsNetwork,
    replicates_field = 25, replicates_genetic = 25)
    }
}

\examples{
data(CmonsRare)
str(CmonsRare) # Structure of data
}

\keyword{ datasets }
