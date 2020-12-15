\name{dhalfnormal}
\alias{dhalfnormal}
\alias{phalfnormal}
\alias{qhalfnormal}
\alias{rhalfnormal}
\alias{ehalfnormal}
\alias{vhalfnormal}
\alias{dhalft}
\alias{phalft}
\alias{qhalft}
\alias{rhalft}
\alias{ehalft}
\alias{vhalft}
\alias{dhalfcauchy}
\alias{phalfcauchy}
\alias{qhalfcauchy}
\alias{rhalfcauchy}
\alias{ehalfcauchy}
\alias{vhalfcauchy}
\title{Half-normal, half-Student-t and half-Cauchy distributions.
}
\description{
  Half-normal, half-Student-t and half-Cauchy density,
  distribution, quantile functions, random number generation,
  and expectation and variance.
}
\usage{
  dhalfnormal(x, scale=1, log=FALSE)
  phalfnormal(q, scale=1)
  qhalfnormal(p, scale=1)
  rhalfnormal(n, scale=1)
  ehalfnormal(scale=1)
  vhalfnormal(scale=1)

  dhalft(x, df, scale=1, log=FALSE)
  phalft(q, df, scale=1)
  qhalft(p, df, scale=1)
  rhalft(n, df, scale=1)
  ehalft(df, scale=1)
  vhalft(df, scale=1)

  dhalfcauchy(x, scale=1, log=FALSE)
  phalfcauchy(q, scale=1)
  qhalfcauchy(p, scale=1)
  rhalfcauchy(n, scale=1)
  ehalfcauchy(scale=1)
  vhalfcauchy(scale=1)
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
  and \sQuote{\code{rhalfnormal()}} generates random deviates.
  The \sQuote{\code{ehalfnormal()}} and \sQuote{\code{vhalfnormal()}}
  functions return the corresponding half-normal distribution's
  expectation and variance, respectively.
  For the
  \sQuote{\code{dhalft()}}, \sQuote{\code{dhalfcauchy()}} and related
  function it works analogously.
}
\references{
  C. Roever, R. Bender, S. Dias, C.H. Schmid, H. Schmidli, S. Sturtz,
  S. Weber, T. Friede.
  On weakly informative prior distributions for the heterogeneity
  parameter in Bayesian random-effects meta-analysis.
  \emph{\href{https://arxiv.org/abs/2007.08352}{arXiv preprint 2007.08352}}
  (submitted for publication), 2020.
  
  A. Gelman. Prior distributions for variance parameters in hierarchical
  models. \emph{Bayesian Analysis}, \bold{1}(3):515-534, 2006.

  F. C. Leone, L. S. Nelson, R. B. Nottingham.
  \href{https://doi.org/10.2307/1266560}{The folded normal
  distribution}. \emph{Technometrics}, \bold{3}(4):543-550, 1961. 

  N. G. Polson, J. G. Scott.
  On the half-Cauchy prior for a global scale parameter.
  \emph{Bayesian Analysis}, \bold{7}(4):887-902, 2012.

  S. Psarakis, J. Panaretos.
  \href{https://doi.org/10.1080/03610929008830342}{The folded t
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
