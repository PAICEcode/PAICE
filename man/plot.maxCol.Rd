\name{plot.maxCol}

\alias{plot.maxCol}

\title{Plot asymptotic estimators of colonization events}

\description{
    Plots for the estimators calculated by \code{\link{maxCol}}.
}

\usage{
\method{plot}{maxCol}(x, xlim, ylim, col, xlabbotton, xlabtop, ylab, main,
     pch = 16, lty = 1, lwd = 2, cex = 1, estimation = TRUE,
     legend = TRUE, \dots)
}

\arguments{
    \item{x}{
        an object of \code{\link[base]{class}} \code{"maxCol"} returned by
        \code{\link{maxCol}} function.
    }
    \item{xlim, ylim}{
        numeric vector containing limits of x and y axis of the plot (min, max).
    }
    \item{col}{
        character vector containing colour of both estimation: genetic and field.
    }
    \item{xlabbotton}{
        a title of the x axis at the bottom of the plot. It correspond with
        genetic estimation.
    }
    \item{xlabtop}{
        a title of the x axis at the top of the plot. It correponds with field
        estimation.
    }
    \item{ylab}{
        a title of the y axis of the plot.
    }
    \item{main}{
        an overall title for the plot.
    }
    \item{pch}{
        indicate symbol used for points of the plot (by default
        \code{pch = 16}). See \code{\link[graphics]{par}} for additional
        information.
    }
    \item{lty}{
        type of lines used in the plot for curve fitting (by default
        \code{lty = 1}). See \code{\link[graphics]{par}} for additional
        information.
    }
    \item{lwd}{
        width of lines used in the plot for curve fitting (by default
        \code{lwd = 2}). See \code{\link[graphics]{par}} for additional
        information.
    }
    \item{cex}{
        size of elements in the plot.
    }
    \item{estimation}{
        logical. If it is \code{TRUE} an estimation of estimation of number of
        colonization events is plotted at the right size of the plot.
    }
    \item{legend}{
        logial. If it is \code{TRUE} the legend is added.
    }
    \item{\dots}{
        aditional graphical parameters (see \code{\link[graphics]{par}}) for
        aditional information.
    }
}

\details{
    Genetic and field estimation are fitted to Michaelis-Menten equation
    following Coello et al. (2022).
}

\value{
    Plot returned by this function represent estimations calculated by
    \code{\link{maxCol}}. The two curves representing both estimators: genetic
    and field. Each point represent the mean of number of colonization events
    inferred by all replicates at this sampling level. Curves represent
    Michaelis-Menten equation fitted to this dataset. If it is plotted, right
    side of the plot represent the number of colonization events estimated by
    this fitting curve for each estimation, including the conficende interval of
    this estimation.
}

\references{
    Coello, A.J., Fernandez-Mazuecos, M., Heleno, R.H., Vargas, P. (2022).
    PAICE: A new R package to estimate the number of inter-island colonizations
    considering haplotype data and sample size. \emph{Journal of Biogeography},
    49(4), 577-589.DOI:
    \href{https://onlinelibrary.wiley.com/doi/10.1111/jbi.14341}{10.1111/jbi.14341}
}

\seealso{
    \code{\link{maxCol}} to fit the accumulation curve of colonization events
    and estimate the number of colonization events.
}

\examples{
# Use 'CmonsRare' data, a dataset generated using 25 replicates
# in both genetic and field sampling
data(CmonsRare)
maxcol <- maxCol(data = CmonsRare)
plot(maxcol)
}
\keyword{ aplot }
