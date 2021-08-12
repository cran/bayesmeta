\name{convolve}
\alias{convolve}
\title{Convolution of two probability distributions}
\description{
  Compute the convolution of two probability distributions, specified
  through their densities or cumulative distribution functions (CDFs).
}
\usage{
  convolve(dens1, dens2,
           cdf1=Vectorize(function(x){integrate(dens1,-Inf,x)$value}),
           cdf2=Vectorize(function(x){integrate(dens2,-Inf,x)$value}),
           delta=0.01, epsilon=0.0001)
}
\arguments{
  \item{dens1, dens2}{the two distributions' probability density functions.}
  \item{cdf1, cdf2}{the two distributions' cumulative distribution functions.}
  \item{delta, epsilon}{the parameters specifying the desired accuracy
    for approximation of the convolution, and with that determining the
    number of support points being used internally. Smaller values
    imply greater accuracy and greater computational burden (Roever and
    Friede, 2017).
}
}
\details{
  The distribution of the \emph{sum} of two (independent) random
  variables technically results as a \emph{convolution} of their
  probability distributions. In some cases, the calculation of
  convolutions may be done analytically; e.g., the sum of two normally
  distributed random variables again turns out as normally distributed
  (with mean and variance resulting as the sums of the original ones).
  In other cases, convolutions may need to be determined
  numerically. One way to achieve this is via the \emph{DIRECT}
  algorithm; the present implementation is the one discussed by Roever
  and Friede (2017). Accuracy of the computations is determined by the
  \code{delta} (maximum divergence \eqn{\delta}) and \code{epsilon}
  (tail probability \eqn{\epsilon}) parameters.

  Convolutions here are used within the \code{\link{funnel}()} function (to
  generate funnel plots), but are often useful more generally. The
  original probability distributions may be specified via their
  probability density functions or their cumulative distribution
  functions (CDFs). The \code{convolve()} function returns the
  convolution's density, CDF and quantile function (inverse CDF).
}
\value{
  A \code{list} with elements
  \item{delta}{the \eqn{\delta} parameter.}
  \item{epsilon}{the \eqn{\epsilon} parameter.}
  \item{binwidth}{the bin width.}
  \item{bins}{the total number of bins.}
  \item{support}{a \code{matrix} containing the support points used
    internally for the convolution approximation.}
  \item{density}{the probability density function.}
  \item{cdf}{the cumulative distribution function (CDF).}
  \item{quantile}{the quantile function (inverse CDF).}
}
\references{
  C. Roever, T. Friede.
  Discrete approximation of a mixture distribution via restricted divergence.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{26}(1):217-222, 2017.
  \doi{10.1080/10618600.2016.1276840}.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{bayesmeta}}, \code{\link{funnel}}.
}
\examples{
\dontrun{
#  Skew-normal / logistic example:

dens1 <- function(x, shape=4)
# skew-normal distribution's density
# see also: http://azzalini.stat.unipd.it/SN/Intro
{
  return(2 * dnorm(x) * pnorm(shape * x))
}

dens2 <- function(x)
# logistic distribution's density
{
  return(dlogis(x, location=0, scale=1))
}

rskewnorm <- function(n, shape=4)
# skew-normal random number generation
# (according to http://azzalini.stat.unipd.it/SN/faq-r.html)
{
  delta <- shape / sqrt(shape^2+1)
  u1 <- rnorm(n); v <- rnorm(n)
  u2 <- delta * u1 + sqrt(1-delta^2) * v
  return(apply(cbind(u1,u2), 1, function(x){ifelse(x[1]>=0, x[2], -x[2])}))
}

# compute convolution:
conv <- convolve(dens1, dens2)

# illustrate convolution:
n <- 100000
x <- rskewnorm(n)
y <- rlogis(n)
z <- x + y

# determine empirical and theoretical quantiles:
p      <- c(0.001,0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
equant <- quantile(z, prob=p)
tquant <- conv$quantile(p)

# show numbers:
print(cbind("p"=p, "empirical"=equant, "theoretical"=tquant))

# draw Q-Q plot:
rg <- range(c(equant, tquant))
plot(rg, rg, type="n", asp=1, main="Q-Q-plot",
     xlab="theoretical quantile", ylab="empirical quantile")
abline(0, 1, col="grey")
points(tquant, equant)
}
}
\keyword{ math }
