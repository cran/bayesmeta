\name{dinvchi}
\alias{dinvchi}
\alias{pinvchi}
\alias{qinvchi}
\alias{rinvchi}
\alias{einvchi}
\alias{vinvchi}
\title{Inverse-Chi distribution.
}
\description{
  (Scaled) inverse-Chi density, distribution, and quantile functions,
  random number generation and expectation and variance.
}
\usage{
  dinvchi(x, df, scale=1, log=FALSE)
  pinvchi(q, df, scale=1, lower.tail=TRUE, log.p=FALSE)
  qinvchi(p, df, scale=1, lower.tail=TRUE, log.p=FALSE)
  rinvchi(n, df, scale=1)
  einvchi(df, scale=1)
  vinvchi(df, scale=1)
}
\arguments{
  \item{x, q}{quantile.}
  \item{p}{probability.}
  \item{n}{number of observations.}
  \item{df}{degrees-of-freedom parameter (\eqn{>0}).}
  \item{scale}{scale parameter (\eqn{>0}).}
  \item{log}{logical; if \code{TRUE}, logarithmic density will be
    returned.}
  \item{lower.tail}{logical; if \code{TRUE} (default),
    probabilities are P(X <= x), otherwise, P(X > x).}
  \item{log.p}{logical; if \code{TRUE}, probabilities p are
    returned as log(p).}
}
\details{
  The \bold{(scaled) inverse-Chi distribution} is defined as the
  distribution of the (scaled) inverse of the square root of a
  Chi-square-distributed random variable. It is a special case of the
  \emph{square-root inverted-gamma} distribution (with
  \eqn{\alpha=\nu/2} and \eqn{\beta=1/2}) (Bernardo and Smith;
  1994). Its probability density function is given by
  \deqn{p(x) \;=\; \frac{2^{(1-\nu/2)}}{s \, \Gamma(\nu/2)}
    \Bigl(\frac{s}{x}\Bigr)^{(\nu+1)}
    \exp\Bigl(-\frac{s^2}{2\,x^2}\Bigr)}{p(x)=2^(1-nu/2)/(s*Gamma(nu/2))
    * (s/x)^(nu+1) * exp(-s^2/(2x^2))}
  where \eqn{\nu} is the \emph{degrees-of-freedom} and \eqn{s} the
  \emph{scale} parameter.
}
\value{
  \sQuote{\code{dinvchi()}} gives the density function,
  \sQuote{\code{pinvchi()}} gives the cumulative distribution
  function (CDF),
  \sQuote{\code{qinvchi()}} gives the quantile function (inverse CDF),
  and \sQuote{\code{rinvchi()}} generates random deviates.
  The \sQuote{\code{einvchi()}} and \sQuote{\code{vinvchi()}}
  functions return the corresponding distribution's
  expectation and variance, respectively.
}
\references{
  J.M. Bernardo, A.F.M. Smith. \emph{Bayesian theory},
  Appendix A.1. Wiley, Chichester, UK, 1994.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{dhalfnormal}},
  \code{\link{dhalft}}.
}
\examples{
#################################
# illustrate Chi^2 - connection;
# generate Chi^2-draws:
chi2 <- rchisq(1000, df=10)
# transform:
invchi <- sqrt(1 / chi2)
# show histogram:
hist(invchi, probability=TRUE, col="grey")
# show density for comparison:
x <- seq(0, 1, length=100)
lines(x, dinvchi(x, df=10, scale=1), col="red")
# compare theoretical and empirical moments:
rbind("theoretical" = c("mean" = einvchi(df=10, scale=1),
                        "var"  = vinvchi(df=10, scale=1)),
      "Monte Carlo" = c("mean" = mean(invchi),
                        "var"  = var(invchi)))

##############################################################
# illustrate the normal/Student-t - scale mixture connection;
# specify degrees-of-freedom:
df <- 5
# generate standard normal draws:
z <- rnorm(1000)
# generate random scalings:
sigma <- rinvchi(1000, df=df, scale=sqrt(df))
# multiply to yield Student-t draws:
t <- z * sigma
# check Student-t distribution via a Q-Q-plot:
qqplot(qt(ppoints(length(t)), df=df), t)
abline(0, 1, col="red")
}
\keyword{ distribution }
