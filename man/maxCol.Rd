\name{maxCol}

\alias{maxCol}

\title{Asymptotic estimation of the number of colonization events}

\description{
    This function calculates the asymptotic estimator of colonization events
    from both curves generated using the {\link{rarecol}} function.
}

\usage{
maxCol(data, level = 0.95, del = 0.05, method = 1)
}

\arguments{
    \item{data}{
        an object of \code{\link[base]{class}} \code{"rarecol"} that contains
        output from \code{\link{rarecol}}.
    }
    \item{level}{
        numeric. Determines the confidence interval used to estimate error in
        Michaelis-Menten equation parameters. By default 0.95.
    }
    \item{del}{
        numeric. Determines the interval of values to be deleted to avoid the
        influence of extreme values. By default 0.05 (i.e. deleted values below
        2.5 quantile and above 97.5 quantile).
    }
    \item{method}{
        numeric. Indicates if the algorithm should try to fit the curve by
        assigning a value to the intercept in genetic rarefaction curves
        (\code{method = 1}) or discard these cases when it is not possible to
        fit the curve with all values (\code{method = 0}). By default,
        \code{method = 1}.
    }
}

\details{
    This function calculates the number of colonization events estimated by both
    resampling methods used in the function \code{\link{rarecol}}. The first
    estimation (genetic estimation) corresponds to resampling first at genetic
    level (number of variable positions) and then, per each variable position, a
    complete resampling of the number of populations is done. The second
    estimation (field estimation) corresponds to the opposite resampling, it is
    done first at field level (number of populations) and then, per each
    population, a complete resampling of the number of variable positions is
    done.

    For each curve, the function first estimates the asymptote (estimated number
    of colonization events) for each level of the second resampling (populations
    in the first estimation and variable positions in the second estimation)
    using the mean value of all replicates at each point. Then, these
    estimations are used to build the final curve estimating the number of
    colonization events for each resampling methodology. This final curve uses
    estimations calculated previously, and the asymptote of the curve is
    calculated by using mean points for each value of the first resampling
    method (variable positions in the first estimation and populations in the
    second estimation). The asymptote is calculated by fitting the curve to a
    Michaelis-Menten equation following Coello et al. (in prep.).

    The confidence interval for the estimated number of colonization events is
    calculated with the \code{\link[stats]{confint}} function. Curve fitting is done using the \code{\link[stats]{nls}} function.
}

\value{
    This function returns an object of \code{\link[base]{class}} \code{"maxCol"}
    consisting in a list of the following elements:

    \item{DataGen}{
        a data frame containing the mean estimated number of colonization events
        per number of variable positions in the genetic estimation.
    }

    \item{FormulaGen}{
        formula used to fit final curve in the genetic estimation.
    }
    \item{DataField}{
        a data frame containing the mean estimated number of colonization events
        per population in the field estimation.
    }
    \item{FormulaField}{
        formula used to fit the final curve in the field estimation.
    }
    \item{Summary}{
        a matrix containing the estimated number of colonization events of each
        estimation (genetic and field). Minimum and maximum are calculated using
        the confidence interval indicated.
    }
    \item{ParametersGen}{
        a matrix containing the value of each parameter to fit a
        Michaelis-Menten equation for genetic estimation. The minimum and
        maximum of each parameter according to the confidence interval
        indicated. This equation is described as: colonization events = M *
        positions / (K + positions) + c.
    }
    \item{ParametersField}{
        a matrix containing the value of each parameter to fit a
        Michaelis-Menten equation for field estimation. The minimum and maximum
        of each parameter according to the confidence interval indicated. This
        equation is described as: colonization events = M * (populations - 1) /
        (K + populations - 1).
    }
    \item{ConfintLevel}{
        a vector containing the confidence interval used to calculate minimum
        and maximum for each parameter.
    }
    \item{DeletedData}{
        a vector containing the interval of extreme values deleted to do the fit
        of the second accumulation curve.
    }

    The function \code{print} returns the number of colonization events inferred
    for each estimation (genetic and field) and the interval of confidence of
    these estimations.The function \code{summary} shows a detailed description
    of parameters used to fit both curves, the formula used to fit these curves
    and the confidence of interval of each parameter.
}

\references{
    Coello, A.J, Heleno, R.H., Fernandez-Mazuecos, M. Vargas, P.
    (in prep.).
    PAICE: a novel method for Phylogeographic Analysis of Island
    Colonization Events considering sampling effort.
}

\note{
    To show a detailed description of inferred colonization events in the most
    complete case use the function \code{\link{colonization}}.
}

\seealso{
    \code{\link{rarecol}} to build rarefaction curves of colonization events. To
    describe the number of colonization events inferred in the most complete
    case use the function \code{\link{colonization}}. \code{\link{plot.maxCol}}
    to plot the result of this function.
}
\examples{
# Use 'CmonsRare' data, a dataset generated using 25 replicates
# in both genetic and field sampling
data(CmonsRare)
maxcol <- maxCol(data = CmonsRare)
maxcol # Number of colonization estimated in each curve
summary(maxcol) # Description of curves
plot(maxcol) # Plotting estimations
# Plot all the information
par(mfrow = c(2, 2))
plot(CmonsRare) # First two plots with rarefaction curves
par(fig = c(0, 1, 0, 0.5), new = TRUE)
plot(maxcol) # Third plot with estimations
}
\keyword{ univar }
