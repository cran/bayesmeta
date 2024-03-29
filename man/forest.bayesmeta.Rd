\name{forest.bayesmeta}
\alias{forest.bayesmeta}
\title{
  Generate a forest plot for a \code{\link{bayesmeta}} object
  (based on the \code{metafor} package's plotting functions).
}
\description{
  Generates a forest plot, showing individual estimates along with their
  95 percent confidence intervals, resulting effect estimate and
  prediction interval.
}
\usage{
  \method{forest}{bayesmeta}(x, xlab="effect size", refline=0, cex=1,...)
}
\arguments{
  \item{x}{
    a \code{\link{bayesmeta}} object.
  }
  \item{xlab}{
    title for the x-axis.
  }
  \item{refline}{
    value at which a vertical \sQuote{reference} line should be drawn
    (default is 0). The line can be suppressed by setting this argument
    to \sQuote{NA}.
  }
  \item{cex}{
    character and symbol expansion factor.
  }
  \item{\dots}{other arguments.
  }
}
\details{
  Generates a simple forest plot illustrating the underlying data and
  resulting estimates (effect estimate and prediction interval).
}
\note{This function requires the \pkg{metafor} package to be installed.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\references{
  C. Lewis and M. Clarke.
  Forest plots: trying to see the wood and the trees.
  \emph{BMJ}, \bold{322}:1479, 2001.
  \doi{10.1136/bmj.322.7300.1479}.

  R.D. Riley, J.P. Higgins and J.J. Deeks.
  Interpretation of random effects meta-analyses.
  \emph{BMJ}, \bold{342}:d549, 2011.
  \doi{10.1136/bmj.d549}.
} 
\seealso{
  \code{\link{bayesmeta}}, \code{\link[metafor]{forest.default}},
  \code{\link[metafor]{addpoly}}, \code{\link{forestplot.bayesmeta}}
}
\examples{
data("CrinsEtAl2014")

\dontrun{
# compute effect sizes (log odds ratios) from count data
# (using "metafor" package's "escalc()" function):
require("metafor")
es.crins <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)
# derive a prior distribution for the heterogeneity:
tp.crins <- TurnerEtAlPrior("surgical", "pharma", "placebo / control")
# perform meta-analysis:
ma.crins <- bayesmeta(es.crins, tau.prior=tp.crins$dprior)

########
# plot:
forest(ma.crins, xlab="log odds ratio")

forest(ma.crins, trans=exp, refline=1, xlab="odds ratio")
}
}
\keyword{ hplot }
