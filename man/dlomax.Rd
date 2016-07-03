\name{dlomax}
\alias{dlomax}
\alias{dhalfnormal}
\alias{dhalft}
\alias{dhalfcauchy}
\alias{dinvchisq}
\alias{drayleigh}
\title{Lomax, half-normal, half-t, half-Cauchy, scaled inverse
  chi-squared and Rayleigh densities.}
\description{
  Compute probability density functions for Lomax, half-normal,
  half-t, half-Cauchy, scaled inverse chi-squared
  and Rayleigh probability distributions.
}
\usage{
  dlomax(x, scale=1, shape=1, log=FALSE)
  dhalfnormal(x, scale=1, log=FALSE)
  dhalft(x, scale=1, df, log=FALSE)
  dhalfcauchy(x, scale=1, log=FALSE)
  dinvchisq(x, scale=1, df=1, log=FALSE)
  drayleigh(x, scale=1, log=FALSE)
}
\arguments{
  \item{x}{argument.}
  \item{shape}{shape parameter (\eqn{>0}).}
  \item{scale}{scale parameter (\eqn{>0}).}
  \item{df}{degrees-of-freedom parameter (\eqn{>0}).}
  \item{log}{logical; if \code{TRUE}, logarithmic density will be returned.}
}
\details{
  The \bold{Lomax distribution} is a special case of a Pareto distribution (2nd
  kind), it is a heavy-tailed distribution that also arises as a mixture
  of an exponential distribution whose  parameter follows a gamma
  distribution. The probability density function of a Lomax distributed
  variable with shape \eqn{\alpha>0} and scale \eqn{\lambda>0} is given by
  \deqn{p(x) = (\alpha / \lambda) (1 + x / \lambda)^{-(\alpha+1)}.}
  Its mean is \eqn{\lambda / (\alpha-1)} (for \eqn{\alpha>1}) and its
  median is \eqn{\alpha(2^{1/\alpha}-1)}{alpha*(2^(1/alpha)-1)}. Its variance is finite only
  for \eqn{\alpha > 2} and equals
  \eqn{(\lambda^2 \alpha) / ((\alpha-1)^2 (\alpha-2))}{(lambda^2*alpha) / ((alpha-1)^2 * (alpha-2))}.

  The \bold{half-normal distribution} is simply a zero-mean normal distribution
  that is restricted to take only positive values. The \emph{scale}
  parameter \eqn{\sigma} here corresponds to the underlying normal
  distribution's standard deviation. Its mean is \eqn{\sigma
  \sqrt{2/\pi}}{sigma*sqrt(2*pi)}, and its variance is \eqn{\sigma^2
  (1-2/\pi)}{sigma^2*(1-2*pi)}. 
  Analogously, the \bold{half-t distribution} is a truncated Student-t
  distribution with \code{df} degrees-of-freedom,
  and the \bold{half-Cauchy distribution} is again a special case of the
  half-t distribution with \code{df=1} degrees of freedom.

  The \bold{scaled inverse \eqn{\chi^2} distribution} constitutes the
  (conjugate) posterior distribution in a simple normal model with
  unknown variance. It posesses two parameters, the \emph{degrees of
  freedom} \eqn{\nu>0} and a \emph{scale} \eqn{s>0}. Its expectation is
  \eqn{(\nu / (\nu+2)) s^2} and the variance is
  \eqn{((2\nu^2)/((\nu-2)^2 (\nu-4))) s^4}.
  \bold{Note} that this distribution is refers to the \emph{variance},
  so when you want to use it as a prior for \eqn{\tau}, a \emph{standard
  deviation} parameter, you need to apply a correction for the square
  root transformation. A scaled inverse \eqn{\chi^2} prior for
  \eqn{\tau^2} with density \eqn{f(x)} corresponds to a prior with
  density \eqn{2\tau f(\tau^2)}{2*tau*f(tau^2)} for \eqn{\tau}.

  The \bold{Rayleigh distribution} arises as the distribution of the
  square root of an exponentially distributed (or
  \eqn{\chi^2_2}-distributed) random variable.
  If \eqn{X} follows an exponential distribution with rate \eqn{\lambda}
  and expectation \eqn{1/\lambda}, then \eqn{Y=\sqrt{X}}{Y=sqrt(X)} follows a
  Rayleigh distribution with scale
  \eqn{\sigma=1/\sqrt{2\lambda}}{sigma=1/sqrt(2*lambda)} and
  expectation \eqn{\sqrt{\pi/(4\lambda)}}{sqrt(pi/(4*lambda))}.
}
\value{
  a vector of probability density values corresponding to \code{x}.
}
\references{
  N.L. Johnson, S. Kotz, N. Balakrishnan. \emph{Continuous univariate
  distributions}, volume 1. Wiley, New York, 2nd edition, 1994.

  A. Gelman, J.B. Carlin, H.S. Stern, D.B. Rubin. \emph{Bayesian data
  analysis}. Chapman & Hall / CRC, Boca Raton, 1997.

  A. Gelman. Prior distributions for variance parameters in hierarchical
  models. \emph{Bayesian Analysis}, \bold{1}(3):515-534, 2006.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{dexp}}, \code{\link{dnorm}}, \code{\link{bayesmeta}}.
}
\examples{
#######################
# illustrate densities:
x <- seq(0,6,le=200)
plot(x, dexp(x, rate=1), type="l", col="orange", ylim=c(0,1),
     xlab=expression(tau), ylab=expression("probability density "*f(tau)))
lines(x, dlomax(x), col="red")
lines(x, dhalfnormal(x), col="blue")
lines(x, dhalft(x, df=3), col="blue", lty="dashed")
lines(x, dhalfcauchy(x), col="blue", lty="dotted")
lines(x, 2*x*dinvchisq(x^2), col="green")
lines(x, drayleigh(x, scale=1/sqrt(2)), col="cyan3")
abline(h=0, v=0, col="grey")

# show log-densities (note the differing tail behaviour):
plot(x, dexp(x, rate=1), type="l", col="orange", ylim=c(0.001,1), log="y",
     xlab=expression(tau), ylab=expression("probability density "*f(tau)))
lines(x, dlomax(x), col="red")
lines(x, dhalfnormal(x), col="blue")
lines(x, dhalft(x, df=3), col="blue", lty="dashed")
lines(x, dhalfcauchy(x), col="blue", lty="dotted")
lines(x, 2*x*dinvchisq(x^2), col="green")
lines(x, drayleigh(x, scale=1/sqrt(2)), col="cyan3")
abline(v=0, col="grey")

#######################################################
# illustrate proper/improper priors using example data:
data("Rubin1981")
\dontrun{
# analysis using improper uniform prior
# (may take a few seconds to compute!):
schools1 <- bayesmeta(y=Rubin1981[,"effect"], sigma=Rubin1981[,"stderr"],
                      labels=Rubin1981[,"school"])
# analysis using weakly informative half-Cauchy prior
# (may take a few seconds to compute!):
schools2 <- bayesmeta(y=Rubin1981[,"effect"], sigma=Rubin1981[,"stderr"],
                      labels=Rubin1981[,"school"],
                      tau.prior=function(x){return(dhalfcauchy(x, scale=25))})

t <- seq(0,50,le=200)
plot(t, schools2$dposterior(tau=t), col="blue", type="l",
     xlab=expression("heterogeneity "*tau), ylab="probability density")
lines(t, schools1$dposterior(tau=t), col="red")
lines(t, schools2$dprior(tau=t), col="blue", lty="dashed")
lines(t, rep(1/max(t),length(t)), col="red", lty="dashed")
legend(max(t), 0.1, xjust=1, yjust=1,
       legend=c("prior A: improper uniform","posterior A",
                "prior B: proper half-Cauchy","posterior B"),
       col=c("red","red","blue","blue"), lty=c(2,1,2,1))
abline(v=0, h=0, col="grey")
}
}
\keyword{ distribution }
