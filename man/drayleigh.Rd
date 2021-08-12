\name{drayleigh}
\alias{drayleigh}
\alias{prayleigh}
\alias{qrayleigh}
\alias{rrayleigh}
\alias{erayleigh}
\alias{vrayleigh}
\title{The Rayleigh distribution.
}
\description{
  Rayleigh density, distribution, quantile function, random number
  generation, and expectation and variance.
}
\usage{
  drayleigh(x, scale=1, log=FALSE)
  prayleigh(q, scale=1)
  qrayleigh(p, scale=1)
  rrayleigh(n, scale=1)
  erayleigh(scale=1)
  vrayleigh(scale=1)
}
\arguments{
  \item{x, q}{quantile.}
  \item{p}{probability.}
  \item{n}{number of observations.}
  \item{scale}{scale parameter (\eqn{>0}).}
  \item{log}{logical; if \code{TRUE}, logarithmic density will be returned.}
}
\details{
  The Rayleigh distribution arises as the distribution of the
  square root of an exponentially distributed (or
  \eqn{\chi^2_2}-distributed) random variable.
  If \eqn{X} follows an exponential distribution with rate \eqn{\lambda}
  and expectation \eqn{1/\lambda}, then \eqn{Y=\sqrt{X}}{Y=sqrt(X)} follows a
  Rayleigh distribution with scale
  \eqn{\sigma=1/\sqrt{2\lambda}}{sigma=1/sqrt(2*lambda)} and
  expectation \eqn{\sqrt{\pi/(4\lambda)}}{sqrt(pi/(4*lambda))}.

  Note that the exponential distribution is the \emph{maximum entropy
    distribution} among distributions supported on the positive real
  numbers and with a pre-specified expectation; so the Rayleigh
  distribution gives the corresponding distribution of its square root.
}
\value{
  \sQuote{\code{drayleigh()}} gives the density function,
  \sQuote{\code{prayleigh()}} gives the cumulative distribution
  function (CDF),
  \sQuote{\code{qrayleigh()}} gives the quantile function (inverse CDF),
  and \sQuote{\code{rrayleigh()}} generates random deviates.
  The \sQuote{\code{erayleigh()}} and \sQuote{\code{vrayleigh()}}
  functions return the corresponding Rayleigh distribution's
  expectation and variance, respectively.
}
\references{
  C. Roever, R. Bender, S. Dias, C.H. Schmid, H. Schmidli, S. Sturtz,
  S. Weber, T. Friede.
  On weakly informative prior distributions for the heterogeneity
  parameter in Bayesian random-effects meta-analysis.
  \emph{Research Synthesis Methods}, \bold{12}(4):448-474, 2021.
  \doi{10.1002/jrsm.1475}.
  
  N.L. Johnson, S. Kotz, N. Balakrishnan. \emph{Continuous univariate
  distributions}, volume 1. Wiley, New York, 2nd edition, 1994.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{dexp}}, \code{\link{dlomax}},
  \code{\link{dhalfnormal}}, \code{\link{dhalft}}, \code{\link{dhalfcauchy}},
  \code{\link{TurnerEtAlPrior}}, \code{\link{RhodesEtAlPrior}},
  \code{\link{bayesmeta}}.
}
\examples{
########################
# illustrate densities:
x <- seq(0,6,le=200)
plot(x, drayleigh(x, scale=0.5), type="l", col="green",
     xlab=expression(tau), ylab=expression("probability density "*f(tau)))
lines(x, drayleigh(x, scale=1/sqrt(2)), col="red")
lines(x, drayleigh(x, scale=1), col="blue")
abline(h=0, v=0, col="grey")

###############################################
# illustrate exponential / Rayleigh connection
# via a quantile-quantile plot (Q-Q-plot):
N <- 10000
exprate <- 5
plot(sort(sqrt(rexp(N, rate=exprate))),
     qrayleigh(ppoints(N), scale=1/sqrt(2*exprate)))
abline(0, 1, col="red")

###############################################
# illustrate Maximum Entropy distributions
# under similar but different constraints:
mu <- 0.5
tau <- seq(0, 4*mu, le=100)
plot(tau, dexp(tau, rate=1/mu), type="l", col="red", ylim=c(0,1/mu),
     xlab=expression(tau), ylab="probability density")
lines(tau, drayleigh(tau, scale=1/sqrt(2*1/mu^2)), col="blue")
abline(h=0, v=0, col="grey")
abline(v=mu, col="darkgrey"); axis(3, at=mu, label=expression(mu))
# explicate constraints:
legend("topright", pch=15, col=c("red","blue"),
       c(expression("Exponential:  E["*tau*"]"==mu),
         expression("Rayleigh:  E["*tau^2*"]"==mu^2)))
}
\keyword{ distribution }
