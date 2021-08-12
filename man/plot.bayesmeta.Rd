\name{plot.bayesmeta}
\alias{plot.bayesmeta}
\title{
  Generate summary plots for a \code{\link{bayesmeta}} object.
}
\description{
  Generates a forest plot, and joint and marginal posterior density plots
  for the two parameters of the random-effects meta-analysis model.
}
\usage{
  \method{plot}{bayesmeta}(x, main=deparse(substitute(x)),
     which=1:4, prior=FALSE, forest.margin=8,
     mulim=c(NA,NA), taulim=c(NA,NA),
     violin=FALSE, ...) 
}
\arguments{
  \item{x}{
    a \code{\link{bayesmeta}} object.
  }
  \item{main}{
    a \code{character} string giving the main title for the plot(s).
  }
  \item{which}{
    an indicator of which plots to generate.
  }
  \item{prior}{
    an indicator whether to also draw the prior density in marginal
    posterior density plots.
  }
  \item{forest.margin}{
    the width of the margin to the left of the forest plot. This may
    require some manual tweaking so that the study labels fit properly.
  }
  \item{mulim, taulim}{
    (optional) ranges of effect (mu) and heterogeneity (tau) values
    to be used for plotting.
  }
  \item{violin}{
    an indicator whether to draw the forest plot as a \dQuote{violin plot}.
  }
  \item{\dots}{other graphical parameters.
  }
}
\details{
  Depending on the value of the \code{which} argument, one or several
  plots are generated, namely
  \enumerate{
    \item a forest plot, including a 95\% credible interval (diamond) and
      a 95\% prediction interval (rectangle) for the effect \eqn{\mu}. The
      shown intervals for \eqn{\mu} are based on posterior medians and
      shortest credible intervals (from \code{x$summary}).
      If \code{violin=TRUE}, the forest plot is plotted as a
      \dQuote{violin plot}, i.e., via Gaussian densities for the
      estimates \eqn{y_i}{y[i]} (and their associated uncertainties),
      and the posterior densities for the effect \eqn{\mu}, and for the
      predictive distribution.
    \item a plot of the joint posterior density of heterogeneity
      (\eqn{\tau}) and effect (\eqn{\mu}). Red lines trace the
      contours of constant density corresponding to approximate 2D
      credible regions (based on a \eqn{\chi^2}-approximation to the
      logarithmic posterior density) as labelled. The credible
      regions are only an approximation based on a
      \sQuote{well-behaved}, unimodal posterior; contour lines are
      omitted if the posterior mode is not finite. Blue lines show the
      conditional mean effect \eqn{\mu} as a function of the
      heterogeneity \eqn{\tau} (solid line) along with conditional
      95\% confidence bounds (dashed lines). Green lines indicate
      marginal medians and shortest 95\% credible intervals for
      \eqn{\tau} and \eqn{\mu}. 
    \item the marginal posterior probability density of the effect
      \eqn{\mu} with median and shortest 95\% credible interval
      indicated. Depending on the \code{prior} argument, a dashed line
      showing the prior density is added. Note that for improper priors
      the scaling is arbitrary and may be inappropriate for the plot.
    \item the marginal posterior probability density of the heterogeneity
      \eqn{\tau} with median and shortest 95\% credible interval
      indicated. Depending on the \code{prior} argument, a dashed line
      showing the prior density is added. Note that for improper priors
      the scaling is arbitrary and may be inappropriate for the plot.
  }
  The joint posterior density plot (2) especially highlights the dependence
  of the effect estimate on the heterogeneity parameter. In a
  \sQuote{conventional} frequentist meta-analysis, one would commonly first
  estimate the heterogeneity \eqn{\tau}, and then fix this value and
  estimate the effect \eqn{\mu} based on the assumption that the
  heterogeneity estimate was the true value. In the joint density plot,
  this would correspond to considering vertical \dQuote{slices} of the
  parameter space, a slice at \eqn{\tau=0} for the fixed-effects model,
  and a slice a a different \eqn{\tau} value for the random-effects
  model, where the blue lines would then indicate the corresponding
  estimate and confidence interval for \eqn{\mu}.
  
  Note that when using the \code{prior=TRUE} argument, the added line
  may end up be outside the plotted range, especially when using
  improper priors with arbitrary normalisation (consider adding it
  \dQuote{manually} instead).
}
\value{
  Returns the supplied \code{bayesmeta} object (\code{x}).
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\references{
  C. Roever.
  Bayesian random-effects meta-analysis using the bayesmeta R package.
  \emph{Journal of Statistical Software}, \bold{93}(6):1-51, 2020.
  \doi{10.18637/jss.v093.i06}.
  
  C. Guddat, U. Grouven, R. Bender and G. Skipka.
  A note on the graphical presentation of prediction intervals in
  random-effects meta-analyses.
  \emph{Systematic Reviews}, \bold{1}(34), 2012.
  \doi{10.1186/2046-4053-1-34}.

  R.D. Riley, J.P. Higgins and J.J. Deeks.
  Interpretation of random effects meta-analyses.
  \emph{BMJ}, \bold{342}:d549, 2011.
  \doi{10.1136/bmj.d549}.
}
\seealso{
  \code{\link{bayesmeta}}, \code{\link{forestplot.bayesmeta}}
}
\examples{
\dontrun{
# example data by Snedecor and Cochran:
data("SnedecorCochran")

# analyze using a weakly informative prior
# (may take a few seconds to compute!):
ma <- bayesmeta(y=SnedecorCochran[,"mean"], sigma=sqrt(SnedecorCochran[,"var"]),
                label=SnedecorCochran[,"no"],
                mu.prior.mean=50, mu.prior.sd=50,
                tau.prior=function(x){dhalfcauchy(x, scale=10)})

# show some plots:
plot(ma, main="Snedecor/Cochran data", prior=TRUE)
plot(ma, main="Snedecor/Cochran data", which=1, violin=TRUE)
}
}
\keyword{ hplot }
