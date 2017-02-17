\name{dlomax}
\alias{dlomax}
\alias{plomax}
\alias{qlomax}
\alias{rlomax}
\title{The Lomax distribution.}
\description{
  Lomax density, distribution and quantile functions,
  and random number generation.
}
\usage{
  dlomax(x, scale=1, shape=1, log=FALSE)
  plomax(q, scale=1, shape=1)
  qlomax(p, scale=1, shape=1)
  rlomax(n, scale=1, shape=1)
}
\arguments{
  \item{x,q}{quantile.}
  \item{p}{probability.}
  \item{n}{number of observations.}
  \item{scale}{scale parameter (\eqn{\lambda > 0}).}
  \item{shape}{shape parameter (\eqn{\alpha > 0}).}
  \item{log}{logical; if \code{TRUE}, logarithmic density will be returned.}
}
\details{
  The Lomax distribution is a heavy-tailed distribution that also is a
  special case of a \emph{Pareto distribution of the 2nd kind}.
  The probability density function of a Lomax distributed variable with
  shape \eqn{\alpha>0} and scale \eqn{\lambda>0} is given by
  \deqn{p(x) = (\alpha / \lambda) (1 + x / \lambda)^{-(\alpha+1)}.}
  The density function is monotonically decreasing in \eqn{x}. Its mean
  is \eqn{\lambda / (\alpha-1)} (for \eqn{\alpha>1}) and its median is
  \eqn{\alpha(2^{1/\alpha}-1)}{alpha*(2^(1/alpha)-1)}. Its variance is
  finite only for \eqn{\alpha > 2} and equals 
  \eqn{(\lambda^2 \alpha) / ((\alpha-1)^2 (\alpha-2))}{(lambda^2*alpha)
    / ((alpha-1)^2 * (alpha-2))}.
  The cumulative distribution function (CDF) is given by
  \deqn{P(x) = 1-(1+ x / \lambda)^{-\alpha}.}

  The Lomax distribution also arises as a \bold{gamma-exponential
  mixture}. Suppose that \eqn{X} is a draw from an exponential
  distribution whose rate \eqn{\theta} again is drawn from a gamma
  distribution with shape \eqn{a} and scale \eqn{s} (so that
  \eqn{\mathrm{E}[\theta]=as}{E[theta]=as}
  and \eqn{\mathrm{Var}(\theta)=as^2}{Var(theta)=as^2},
  or \eqn{\mathrm{E}[1/\theta]=\frac{1}{s(a+1)}}{E[1/theta]=1/(s*(a+1))}
  and \eqn{\mathrm{Var}(1/\theta)=\frac{1}{s^2(a-1)^2(a-2)}}{Var(1/theta)=1/(s^2*(a-1)^2*(a-2))}). 
  Then the marginal distribution of \eqn{X} is Lomax with scale
  \eqn{1/s} and shape \eqn{a}. Consequently, if the moments of
  \eqn{\theta} are given by \eqn{\mathrm{E}[\theta]=\mu} and
  \eqn{\mathrm{Var}(\theta)=\sigma^2}, then \eqn{X} is Lomax distributed
  with shape
  \eqn{\alpha=\left(\frac{\mu}{\sigma}\right)^2}{alpha=(mu/sigma)^2} and
  scale
  \eqn{\lambda=\frac{\mu}{\sigma^2}=\frac{\alpha}{\mu}}{lambda=mu/sigma^2=alpha/mu}. 
  The gamma-exponential connection is also illustrated in an example below.
}
\value{
  \sQuote{\code{dlomax()}} gives the density function,
  \sQuote{\code{plomax()}} gives the cumulative distribution
  function (CDF),
  \sQuote{\code{qlomax()}} gives the quantile function (inverse CDF),
  and \sQuote{\code{rlomax()}} generates random deviates.
}
\references{
  N.L. Johnson, S. Kotz, N. Balakrishnan. \emph{Continuous univariate
  distributions}, volume 1. Wiley, New York, 2nd edition, 1994.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{dexp}},
  \code{\link{dgamma}},
  \code{\link{dhalfnormal}}, \code{\link{dhalft}}, \code{\link{dhalfcauchy}},
  \code{\link{drayleigh}},
  \code{\link{TurnerEtAlPrior}}, \code{\link{RhodesEtAlPrior}},
  \code{\link{bayesmeta}}.
}
\examples{
#######################
# illustrate densities:
x <- seq(0,6,le=200)
plot(x, dexp(x, rate=1), type="l", col="cyan", ylim=c(0,1),
     xlab=expression(tau), ylab=expression("probability density "*f(tau)))
lines(x, dlomax(x), col="orange")
abline(h=0, v=0, col="grey")

# show log-densities (note the differing tail behaviour):
plot(x, dexp(x, rate=1), type="l", col="cyan", ylim=c(0.001,1), log="y",
     xlab=expression(tau), ylab=expression("probability density "*f(tau)))
lines(x, dlomax(x), col="orange")
abline(v=0, col="grey")

######################################################
# illustrate the gamma-exponential mixture connection;
# specify a number of samples:
N <- 10000
# specify some gamma shape and scale parameters
# (via mixing distribution's moments):
expectation <- 2.0
stdev       <- 1.0
gammashape <- (expectation / stdev)^2
gammascale <- stdev^2 / expectation
print(c("expectation"=expectation, "stdev"=stdev,
        "shape"=gammashape, "scale"=gammascale))
# generate gamma-distributed rates:
lambda <- rgamma(N, shape=gammashape, scale=gammascale)
# generate exponential draws according to gamma-rates:
y <- rexp(N, rate=lambda)
# determine Lomax quantiles accordingly parameterized:
x <- qlomax(ppoints(N), scale=1/gammascale, shape=gammashape)
# compare distributions in a Q-Q-plot:
plot(x, sort(y), log="xy", main="quantile-quantile plot",
     xlab="theoretical quantile", ylab="empirical quantile")
abline(0, 1, col="red")
}
\keyword{ distribution }
