\name{plot.rarecol}

\alias{plot.rarecol}

\title{Plot rarefaction curve of colonization events}

\description{
  Plots for the rarefaction curves produced by \code{\link{rarecol}}.
}
\usage{
\method{plot}{rarecol}(x, xlim1, xlim2, ylim, ylim1, ylim2, palette1, palette2, main1,
    main2, xlab1, xlab2, ylab1, ylab2, las1 = 1, las2 = 1,
    cextText = 0.75, legendbar = TRUE, \dots)
}

\arguments{
    \item{x}{
        an object generated by the \code{\link{rarecol}} function.
    }
    \item{xlim1, xlim2}{
        x limits (min, max) of the two plots.
    }
    \item{ylim1, ylim2}{
        y limits (min, max) of the two plots.
    }
    \item{ylim}{
        y limits (min, max) of the two plots simultaneously. If \code{ylim} is
        defined, the function does not consider \code{ylim1} and \code{ylim2}.
    }
    \item{palette1, palette2}{
        vector of color for lines in plot 1 and plot 2.
    }
    \item{main1, main2}{
        overall title of plot 1 and plot 2.
    }
    \item{xlab1, xlab2}{
        label of x axis of plot 1 and plot 2.
    }
    \item{ylab1, ylab2}{
        label of y axis of plot 1 and plot 2
    }
    \item{las1, las2}{
        numeric. Corresponds to the style of axis labels in plot 1 and plot 2.
        Values: 0 (always parallel to the axis), 1 (always horizontal), 2
        (always perpendicular to the axis), 3 (always vertical). See
        \code{\link[graphics]{par}} for more information.
    }
    \item{cextText}{
        size of legend text.
    }
    \item{legendbar}{
    logical. If \code{TRUE}, it shows a legend bar indicating the color value of
    each variable position.
    }
    \item{\dots}{
        aditional graphical parameters (see \code{\link[graphics]{par}}) for
        additional information.
    }
}

\details{
    The first plot corresponds to the genetic estimation. This plot shows
    accumulation of colonization events as a function of population number. Each
    curve was created for each number of variable positions in the dataset.

    The second plot corresponds to the field estimation. This plot shows
    accumulation of colonization events as a function of the number of variable
    positions. Each curve is created for each number of populations in the
    dataset.
}

\value{
    This function returns two plots corresponding to the two resampling methods
    used in \code{\link{rarecol}}. The first curve corresponds to the "genetic
    estimation" in which a genetic resampling of every possible number of
    variable position is done and, for each resample, a complete resampling of
    population is done. The second curve represents the opposite method
    corresponding to the "field estimation": it first resamples every possible
    number of populations and, for each case, a complete resampling of
    variable positions is done.
}

\references{
    Coello, A.J., Fernandez-Mazuecos, M., Heleno, R.H., Vargas, P. (2022).
    PAICE: A new R package to estimate the number of inter-island colonizations
    considering haplotype data and sample size. \emph{Journal of Biogeography},
    49(4), 577-589.DOI:
    \href{https://onlinelibrary.wiley.com/doi/10.1111/jbi.14341}{10.1111/jbi.14341}
}

\seealso{
    \code{\link{rarecol}} to build a rarefection curve of colonization events.
}
\examples{
# Use 'CmonsRare' data, a dataset generated using 25 replicates
# in both genetic and field sampling
data(CmonsRare)
plot(CmonsRare)
}

\keyword{ aplot }
