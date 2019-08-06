\name{dhalflogistic}
\alias{dhalflogistic}
\alias{phalflogistic}
\alias{qhalflogistic}
\alias{rhalflogistic}
\alias{ehalflogistic}
\alias{vhalflogistic}
\title{Half-logistic distribution.
}
\description{
  Half-logistic density, distribution, and quantile functions,
  random number generation and expectation and variance.
}
\usage{
  dhalflogistic(x, scale=1, log=FALSE)
  phalflogistic(q, scale=1)
  qhalflogistic(p, scale=1)
  rhalflogistic(n, scale=1)
  ehalflogistic(scale=1)
  vhalflogistic(scale=1)
}
\arguments{
  \item{x, q}{quantile.}
  \item{p}{probability.}
  \item{n}{number of observations.}
  \item{scale}{scale parameter (\eqn{>0}).}
  \item{log}{logical; if \code{TRUE}, logarithmic density will be returned.}
}
\details{
  The \bold{half-logistic distribution} is simply a zero-mean logistic distribution
  that is restricted to take only positive values.
  If \eqn{X\sim\mathrm{logistic}}{X ~ logistic}, then
  \eqn{|sX|\sim\mathrm{halflogistic}(\mathrm{scale}\!=\!s)}{|sX| ~ halflogistic(scale=s)}.
}
\value{
  \sQuote{\code{dhalflogistic()}} gives the density function,
  \sQuote{\code{phalflogistic()}} gives the cumulative distribution
  function (CDF),
  \sQuote{\code{qhalflogistic()}} gives the quantile function (inverse CDF),
  and \sQuote{\code{rhalflogistic()}} generates random deviates.
  The \sQuote{\code{ehalflogistic()}} and \sQuote{\code{vhalflogistic()}}
  functions return the corresponding half-logistic distribution's
  expectation and variance, respectively.
}
\references{
  N.L. Johnson, S. Kotz, N. Balakrishnan. \emph{Continuous univariate
  distributions}, volume 2, chapter 23.11. Wiley, New York, 2nd edition, 1994.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{dlogis}}, \code{\link{dhalfnormal}},
  \code{\link{dlomax}}, \code{\link{drayleigh}},
  \code{\link{TurnerEtAlPrior}}, \code{\link{RhodesEtAlPrior}},
  \code{\link{bayesmeta}}.
}
\examples{
#######################
# illustrate densities:
x <- seq(0,6,le=200)
plot(x, dhalfnormal(x), type="l", col="red", ylim=c(0,1),
     xlab=expression(tau), ylab=expression("probability density "*f(tau)))
lines(x, dhalflogistic(x), col="green3")
lines(x, dhalfcauchy(x), col="blue")
lines(x, dexp(x), col="cyan")
abline(h=0, v=0, col="grey")

# show log-densities (note the differing tail behaviour):
plot(x, dhalfnormal(x), type="l", col="red", ylim=c(0.001,1), log="y",
     xlab=expression(tau), ylab=expression("probability density "*f(tau)))
lines(x, dhalflogistic(x), col="green3")
lines(x, dhalfcauchy(x), col="blue")
lines(x, dexp(x), col="cyan")
abline(v=0, col="grey")
}
\keyword{ distribution }
