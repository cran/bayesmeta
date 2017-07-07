\name{dhalfnormal}
\alias{dhalfnormal}
\alias{phalfnormal}
\alias{qhalfnormal}
\alias{rhalfnormal}
\alias{dhalft}
\alias{phalft}
\alias{qhalft}
\alias{rhalft}
\alias{dhalfcauchy}
\alias{phalfcauchy}
\alias{qhalfcauchy}
\alias{rhalfcauchy}
\title{Half-normal, half-Student-t and half-Cauchy distributions.
}
\description{
  Half-normal, half-Student-t and half-Cauchy density,
  distribution, quantile functions and random number generation.
}
\usage{
  dhalfnormal(x, scale=1, log=FALSE)
  phalfnormal(q, scale=1)
  qhalfnormal(p, scale=1)
  rhalfnormal(n, scale=1)

  dhalft(x, scale=1, df, log=FALSE)
  phalft(q, scale=1, df)
  qhalft(p, scale=1, df)
  rhalft(n, scale=1, df)

  dhalfcauchy(x, scale=1, log=FALSE)
  phalfcauchy(q, scale=1)
  qhalfcauchy(p, scale=1)
  rhalfcauchy(n, scale=1)
}
\arguments{
  \item{x, q}{quantile.}
  \item{p}{probability.}
  \item{n}{number of observations.}
  \item{scale}{scale parameter (\eqn{>0}).}
  \item{df}{degrees-of-freedom parameter (\eqn{>0}).}
  \item{log}{logical; if \code{TRUE}, logarithmic density will be returned.}
}
\details{
  The \bold{half-normal distribution} is simply a zero-mean normal distribution
  that is restricted to take only positive values. The \emph{scale}
  parameter \eqn{\sigma} here corresponds to the underlying normal
  distribution's standard deviation:
  if \eqn{X\sim\mathrm{Normal}(0,\sigma^2)}{X ~ Normal(0,sigma)}, then
  \eqn{|X|\sim\mathrm{halfNormal}(\mathrm{scale}\!=\!\sigma)}{|X| ~ halfNormal(scale=sigma)}.
  Its mean is \eqn{\sigma
  \sqrt{2/\pi}}{sigma*sqrt(2/pi)}, and its variance is \eqn{\sigma^2
  (1-2/\pi)}{sigma^2*(1-2/pi)}. 
  Analogously, the \bold{half-t distribution} is a truncated Student-t
  distribution with \code{df} degrees-of-freedom,
  and the \bold{half-Cauchy distribution} is again a special case of the
  half-t distribution with \code{df=1} degrees of freedom.

  Note that (half-) Student-t and Cauchy distributions arise as
  continuous \emph{mixture distributions} of (half-) normal
  distributions. If
  \deqn{Y|\sigma\;\sim\;\mathrm{Normal}(0,\sigma^2)}{Y|sigma  ~  Normal(0,sigma^2)}
  where the standard deviation is
  \eqn{\sigma = \sqrt{k/X}}{sigma=sqrt(k/X)} and \eqn{X} is drawn from a
  \eqn{\chi^2}-distribution with \eqn{k} degrees of freedom, then the
  marginal distribution of \eqn{Y} is Student-t with \eqn{k} degrees of freedom.
}
\value{
  \sQuote{\code{dhalfnormal()}} gives the density function,
  \sQuote{\code{phalfnormal()}} gives the cumulative distribution
  function (CDF),
  \sQuote{\code{qhalfnormal()}} gives the quantile function (inverse CDF),
  and \sQuote{\code{rhalfnormal()}} generates random deviates. For the
  \sQuote{\code{dhalft()}}, \sQuote{\code{dhalfcauchy()}} and related
  function it works analogously.
}
\references{
  A. Gelman. Prior distributions for variance parameters in hierarchical
  models. \emph{Bayesian Analysis}, \bold{1}(3):515-534, 2006.

  F. C. Leone, L. S. Nelson, R. B. Nottingham.
  \href{http://doi.org/10.2307/1266560}{The folded normal
  distribution}. \emph{Technometrics}, \bold{3}(4):543-550, 1961. 

  S. Psarakis, J. Panaretos.
  \href{http://doi.org/10.1080/03610929008830342}{The folded t
  distribution}. \emph{Communications in Statistics - Theory and
  Methods}, \bold{19}(7):2717-2734, 1990.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{dnorm}}, \code{\link{dt}}, \code{\link{dcauchy}},
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
lines(x, dhalft(x, df=3), col="green")
lines(x, dhalfcauchy(x), col="blue")
lines(x, dexp(x), col="cyan")
abline(h=0, v=0, col="grey")

# show log-densities (note the differing tail behaviour):
plot(x, dhalfnormal(x), type="l", col="red", ylim=c(0.001,1), log="y",
     xlab=expression(tau), ylab=expression("probability density "*f(tau)))
lines(x, dhalft(x, df=3), col="green")
lines(x, dhalfcauchy(x), col="blue")
lines(x, dexp(x), col="cyan")
abline(v=0, col="grey")
}
\keyword{ distribution }
