\name{normalmixture}
\alias{normalmixture}
\title{
  Compute normal mixtures
}
\description{
  This function allows to derive density, distribution function and
  quantile function of a normal mixture with fixed mean (\eqn{\mu}) and
  random standard deviation (\eqn{\sigma}).
}
\usage{
  normalmixture(density,
                cdf = Vectorize(function(x){integrate(density,0,x)$value}),
                mu = 0, delta = 0.01, epsilon = 0.0001,
                rel.tol.integrate = 2^16*.Machine$double.eps,
                abs.tol.integrate = rel.tol.integrate,
                tol.uniroot = rel.tol.integrate)
}
\arguments{
  \item{density}{
    the \eqn{\sigma} mixing distribution's probability density
    \code{function}.
  }
  \item{cdf}{
    the \eqn{\sigma} mixing distribution's cumulative distribution
    \code{function}.
  }
  \item{mu}{
    the normal mean (\eqn{\mu}).
  }
  \item{delta, epsilon}{
    the parameters specifying the desired accuracy for approximation of
    the mixing distribution, and with that determining the number of
    \eqn{\sigma} support points being used internally. Smaller values
    imply greater accuracy and greater computational burden (Roever and
    Friede, 2017).
  }
  \item{rel.tol.integrate, abs.tol.integrate, tol.uniroot}{
    the \code{rel.tol}, \code{abs.tol} and \code{tol}
    \sQuote{accuracy} arguments that are passed to
    the \code{\link{integrate}()} or \code{\link{uniroot}()} functions
    for internal numerical integration or root finding
    (see also the help there).
  }
}
\details{
  When a normal random variable
  \deqn{X|\mu,\sigma \;\sim\; \mathrm{Normal}(\mu,\sigma^2)}{X|mu,sigma  ~  Normal(mu, sigma^2)}
  has a fixed mean \eqn{\mu}, but a random standard deviation
  \deqn{\sigma|\phi \;\sim\; \mathrm{G}(\phi)}{sigma|phi  ~  G(phi)}
  following some probability distribution \eqn{\mathrm{G}(\phi)}{G(phi)},
  then the \emph{marginal distribution} of \eqn{X},
  \deqn{X|\mu,\phi}{X|mu,phi}
  is a \emph{mixture distribution} (Lindsay, 1995; Seidel, 2010).

  The mixture distribution's probability density function etc. result
  from integration and often are not available in analytical form.
  The \code{normalmixture()} function approximates density,
  distribution function and quantile function to some pre-set accuracy
  by a \emph{discrete} mixture of normal distributions based on a finite
  number of \eqn{\sigma} values using the \sQuote{DIRECT} algorithm
  (Roever and Friede, 2017).

  Either the \dQuote{\code{density}} or \dQuote{\code{cdf}} argument
  needs to be supplied. If only \dQuote{\code{density}} is given, then
  the CDF is derived via integration, if only \dQuote{\code{cdf}} is
  supplied, then the density function is not necessary.

  In \bold{meta-analysis} applications, mixture distributions arise
  e.g. in the
  context of \bold{prior predictive distributions}. Assuming
  that a study-specific effect \eqn{\theta_i}{theta[i]} \emph{a
    priori} is distributed as
  \deqn{\theta_i|\mu,\tau \;\sim\; \mathrm{Normal}(\mu,\tau^2)}{theta[i]|mu,tau  ~  Normal(mu, tau^2)}
  with a prior distribution for the heterogeneity \eqn{\tau},
  \deqn{\tau|\phi \;\sim\; \mathrm{G}(\phi)}{tau|phi  ~  G(phi)}
  yields a setup completely analogous to the above one.
    
  Since it is sometimes hard to judge what constitutes a sensible
  heterogeneity prior, it is often useful to inspect the implications of
  certain settings in terms of the corresponding \emph{prior predictive
  distribution} of
  \deqn{\theta_i|\mu,\phi}{theta[i]|mu,phi}
  indicating the \emph{a priori} implied variation between studies due
  to heterogeneity alone based on a certain prior setup (Spiegelhalter
  et al., 2004, Sec. 5.7.3). Some examples using different heterogeneity
  priors are illustrated below.
}
\value{
  A \code{list} containing the following elements:
  \item{delta, epsilon}{The supplied design parameters.}
  \item{mu}{the normal mean.}
  \item{bins}{the number of bins used.}
  \item{support}{the matrix containing lower, upper and reference points
    for each bin and its associated probability.}
  \item{density}{the mixture's density \code{function(x)}.}
  \item{cdf}{the mixture's cumulative distribution \code{function(x)} (CDF).}
  \item{quantile}{the mixture's quantile \code{function(p)} (inverse CDF).}
  \item{mixing.density}{the mixing distribution's density
    \code{function()} (if supplied).}
  \item{mixing.cdf}{the mixing distribution's cumulative distribution
    \code{function()}.}
}
\references{
  B.G. Lindsay.
  \emph{Mixture models: theory, geometry and applications}.
  Vol. 5 of \emph{CBMS Regional Conference Series in Probability and
    Statistics}, 
  Institute of Mathematical Statistics, Hayward, CA, USA, 1995.

  C. Roever, T. Friede.
  \href{http://doi.org/10.1080/10618600.2016.1276840}{Discrete
  approximation of a mixture distribution via restricted divergence}.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{26}(1):217-222, 2017.

  C. Roever.
  \href{http://www.doi.org/10.18637/jss.v093.i06}{Bayesian random-effects meta-analysis using the bayesmeta R package}.
  \emph{Journal of Statistical Software}, \bold{93}(6):1-51, 2020.
  
  W.E. Seidel. Mixture models.
  In M. Lovric (ed.) \emph{International Encyclopedia of Statistical
  Science}, Springer, Heidelberg, pp. 827-829, 2010.

  D.J. Spiegelhalter, K.R. Abrams, J.P.Myles. \emph{Bayesian approaches
  to clinical trials and health-care evaluation}. Wiley & Sons, 2004.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{bayesmeta}}.
}
\examples{
##################################################################
# compare half-normal mixing distributions with different scales:
nm05 <- normalmixture(cdf=function(x){phalfnormal(x, scale=0.5)})
nm10 <- normalmixture(cdf=function(x){phalfnormal(x, scale=1.0)})
# (this corresponds to the case of assuming a half-normal prior
# for the heterogeneity tau)

# check the structure of the returned object:
str(nm05)

# show density functions:
# (these would be the marginal (prior predictive) distributions
# of study-specific effects theta[i])
x <- seq(-1, 3, by=0.01)
plot(x, nm05$density(x), type="l", col="blue", ylab="density")
lines(x, nm10$density(x), col="red")
abline(h=0, v=0, col="grey")

# show cumulative distributions:
plot(x, nm05$cdf(x), type="l", col="blue", ylab="CDF")
lines(x, nm10$cdf(x), col="red")
abline(h=0:1, v=0, col="grey")

# determine 5 percent and 95 percent quantiles:
rbind("HN(0.5)"=nm05$quantile(c(0.05,0.95)),
      "HN(1.0)"=nm10$quantile(c(0.05,0.95)))


##################################################################
# compare different mixing distributions
# (half-normal, half-Cauchy, exponential and Lomax):
nmHN <- normalmixture(cdf=function(x){phalfnormal(x, scale=0.5)})
nmHC <- normalmixture(cdf=function(x){phalfcauchy(x, scale=0.5)})
nmE  <- normalmixture(cdf=function(x){pexp(x, rate=2)})
nmL  <- normalmixture(cdf=function(x){plomax(x, shape=4, scale=2)})

# show densities (logarithmic y-axis):
x <- seq(-1, 3, by=0.01)
plot(x,  nmHN$density(x), col="green",  type="l", ylab="density", ylim=c(0.005, 6.5), log="y")
lines(x, nmHC$density(x), col="red")
lines(x, nmE$density(x),  col="blue")
lines(x, nmL$density(x),  col="cyan")
abline(v=0, col="grey")

# show CDFs:
plot(x,  nmHN$cdf(x), col="green",  type="l", ylab="CDF", ylim=c(0,1))
lines(x, nmHC$cdf(x), col="red")
lines(x, nmE$cdf(x),  col="blue")
lines(x, nmL$cdf(x),  col="cyan")
abline(h=0:1, v=0, col="grey")
# add "exponential" x-axis at top:
axis(3, at=log(c(0.5,1,2,5,10,20)), lab=c(0.5,1,2,5,10,20))
# show 95 percent quantiles:
abline(h=0.95, col="grey", lty="dashed")
abline(v=nmHN$quantile(0.95), col="green", lty="dashed")
abline(v=nmHC$quantile(0.95), col="red", lty="dashed")
abline(v=nmE$quantile(0.95),  col="blue", lty="dashed")
abline(v=nmL$quantile(0.95),  col="cyan", lty="dashed")
rbind("half-normal(0.5)"=nmHN$quantile(0.95),
      "half-Cauchy(0.5)"=nmHC$quantile(0.95),
      "exponential(2.0)"=nmE$quantile(0.95),
      "Lomax(4,2)"      =nmL$quantile(0.95))


#####################################################################
# a normal mixture distribution example where the solution
# is actually known analytically: the Student-t distribution.
# If  Y|sigma ~ N(0,sigma^2),  where  sigma = sqrt(k/X)
# and  X|k ~ Chi^2(df=k),
# then the marginal  Y|k  is Student-t with k degrees of freedom.

# define CDF of sigma:
CDF <- function(sigma, df){pchisq(df/sigma^2, df=df, lower.tail=FALSE)}

# numerically approximate normal mixture (with k=5 d.f.):
k <- 5
nmT1 <- normalmixture(cdf=function(x){CDF(x, df=k)})
# in addition also try a more accurate approximation:
nmT2 <- normalmixture(cdf=function(x){CDF(x, df=k)}, delta=0.001, epsilon=0.00001)
# check: how many grid points were required?
nmT1$bins
nmT2$bins

# show true and approximate densities:
x <- seq(-2,10,le=400)
plot(x, dt(x, df=k), type="l")
abline(h=0, v=0, col="grey")
lines(x, nmT1$density(x), col="red", lty="dashed")
lines(x, nmT2$density(x), col="blue", lty="dotted")

# show ratios of true and approximate densities:
plot(x, nmT1$density(x)/dt(x, df=k), col="red",
     type="l", log="y", ylab="density ratio")
abline(h=1, v=0, col="grey")
lines(x, nmT2$density(x)/dt(x, df=k), col="blue")
}
\keyword{ distribution }
