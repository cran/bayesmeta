\name{weightsplot}
\alias{weightsplot}
\alias{weightsplot.default}
\alias{weightsplot.bayesmeta}
\title{
  Illustrate the posterior mean weights for a \code{\link{bayesmeta}} object.
}
\description{
  Generates a bar plot showing individual estimates' posterior mean
  weights, either for the overall mean estimate, or for a shrinkage
  estimate.
}
\usage{
  weightsplot(x, ...)
  \method{weightsplot}{bayesmeta}(x, individual=FALSE, ordered=TRUE,
            extramargin=4, priorlabel="prior mean", main, xlim, ...)
}
\arguments{
  \item{x}{
    a \code{\link{bayesmeta}} object.
  }
  \item{individual}{
    this argument allows to request weights for individual
    \emph{shrinkage estimates}. If \code{FALSE} (the default), weights
    for the overall mean are returned. Otherwise, it may be an integer
    number (\code{1,...,k}) giving the index, or a character string giving
    the label.
  }
  \item{ordered}{
    a logical flag indicating whether to sort weights by their magnitude.
  }
  \item{extramargin}{
    an additional margin to be added to the left side of the plot, in
    order to accomodate the estimates' labels. The value will be added
    to the 2nd element of the margin settings given by
    \sQuote{\code{par("mar")}}. In case study labels still extend beyond
    the figure margin, try increasing this number. See also the
    \code{\link[graphics]{par}()} function's help.
  }
  \item{priorlabel}{
    the label for the effect prior's weight. Only relevant for proper
    effect priors.
  }
  \item{main}{
    the plot's main title.
  }
  \item{xlim}{
    the x-axis range.
  }
  \item{...}{
    other arguments passed on to the
    \code{\link[graphics]{barplot}()} function.
  }
}
\details{
  The individual estimates' contributions to the overall mean estimate
  are commonly illustrated in terms of \emph{weights}, as the resulting
  overall estimate may be expressed as a weighted average of the
  estimates contributing to the analysis. The notion of \dQuote{study
    weights} may also be extended to the Bayesian setting, where these
  result as \emph{posterior mean weights}. Analogous weights may also be
  derived for \emph{shrinkage estimates} (Roever and Friede, 2021).
  
  This function generates a simple bar plot illustrating the
  posterior mean weights. The actual numbers are taken from the
  \code{bayesmeta} object's \dQuote{\code{$weights}} or
  \dQuote{\code{$weights.theta}} elements.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\references{
  C. Roever, T. Friede.
  Bounds for the weight of external data in shrinkage estimation.
  \emph{Biometrical Journal}, \bold{65}(5):1131-1143, 2021.
  \doi{10.1002/bimj.202000227}.
} 
\seealso{
  \code{\link{bayesmeta}}.
}
\examples{
# load example data:
data("Peto1980")
\dontrun{
# compute effect sizes (log odds ratios) from count data:
require("metafor")
peto.es <- escalc(measure="OR",
                  ai=treat.events,   n1i=treat.cases,
                  ci=control.events, n2i=control.cases,
                  slab=publication, data=Peto1980)

# perform meta-analysis:
ma01 <- bayesmeta(peto.es)
# show data and results:
forestplot(ma01)

# check out weights:
ma01$weights
ma01$weights.theta

# illustrate weights:
weightsplot(ma01)
weightsplot(ma01, ordered=FALSE)
weightsplot(ma01, ordered=FALSE, individual="BrMedJ1974")
}
}
\keyword{ hplot }
