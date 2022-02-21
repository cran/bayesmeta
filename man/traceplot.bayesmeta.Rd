\name{traceplot}
\alias{traceplot}
\alias{traceplot.default}
\alias{traceplot.bayesmeta}
\alias{traceplot.bmr}
\title{
  Illustrate conditional means of overall effect as well as study-specific
  estimates as a function of heterogeneity.
}
\description{
  Generates a trace plot of overall mean effect and study-specific
  (shrinkage) estimates as a function of the heterogeneity (\eqn{\tau}). The
  heterogeneity's posterior distribution is also indicated.
}
\usage{
  traceplot(x, ...)
  \method{traceplot}{bayesmeta}(x, mulim, taulim, ci=FALSE,
          rightmargin=8, col=rainbow(x$k), ...)
  \method{traceplot}{bmr}(x, mulim, taulim, ci=FALSE,
          rightmargin=8, col=rainbow(x$k), ...)
}
\arguments{
  \item{x}{
    a \code{\link{bayesmeta}} or \code{\link{bmr}}object.
  }
  \item{mulim, taulim}{(optional) ranges for the effect (mu) and
    heterogeneity (tau) axes. If only one value is given for
    \code{taulim}, then this is taken as the upper limit, and the lower
    limit is zero.
  }
  \item{ci}{
    a logical flag indicating whether to also show (conditional)
    confidence intervals.
  }
  \item{rightmargin}{
    an additional margin to be added to the right side of the plot, in
    order to accomodate the estimates' labels. In case study labels
    still extend beyond the figure margin, try increasing this number.
  }
  \item{col}{
    colors to be used for plotting the (\eqn{k}) estimates.
  }
  \item{...}{
    other arguments passed on to the
    \code{\link[graphics]{plot}()} function.
  }
}
\details{
  For a given heterogeneity (\eqn{\tau}) value, the \emph{conditional}
  posterior distributions of the overall effect (\eqn{\mu}) as well as
  the study-specific parameters (\eqn{\theta_i}{theta[i]}) are again
  normal. The conditional normal moments (mean and variance) then vary
  as functions of the heterogeneity; for large heterogeneity, the
  shrinkage estimates approach the original data (\eqn{\mu_i}{mu[i]}),
  while the overall mean approaches an un-weighted overall average. For
  small heterogeneity, both overall mean as well as study-specific
  estimates are increasingly \emph{shrunk} towards the
  inverse-variance-weighted \sQuote{common-effect} estimate (Roever,
  2020).

  This trace plot illustrates the conditional (overall and
  study-specific) estimates along with the heterogeneity's posterior
  distribution (density) in a layout similar to that utilized
  by Rubin (1981).
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\references{
  C. Roever.
  Bayesian random-effects meta-analysis using the bayesmeta R package.
  \emph{Journal of Statistical Software}, \bold{93}(6):1-51, 2020.
  \doi{10.18637/jss.v093.i06}.
  
  D.B. Rubin.
  Estimation in parallel randomized experiments.
  \emph{Journal of Educational Statistics}, \bold{6}(4):377-401, 1981.
  \doi{10.3102/10769986006004377}.
} 
\seealso{
  \code{\link{bayesmeta}}, \code{\link{bmr}}.
}
\examples{
# load example data:
data("Rubin1981")

\dontrun{
# perform meta-analysis:
ma01 <- bayesmeta(y=Rubin1981[,"effect"], sigma=Rubin1981[,"stderr"],
                  labels=Rubin1981[,"school"], tau.prior="uniform")

# show meta-analysis results:
forestplot(ma01)

# show trace plot:
traceplot(ma01)
}
}
\keyword{ hplot }
