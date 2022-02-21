\name{kldiv}
\alias{kldiv}
\title{Kullback-Leibler divergence of two multivariate normal distributions.}
\description{
  Compute the Kullback-Leiber divergence or
  \emph{symmetrized} KL-divergence based on means
  and covariances of two normal distributions.
}
\usage{
  kldiv(mu1, mu2, sigma1, sigma2, symmetrized=FALSE)
}
\arguments{
  \item{mu1, mu2}{the two mean vectors.}
  \item{sigma1, sigma2}{the two covariance matrices.}
  \item{symmetrized}{logical; if \code{TRUE}, the \emph{symmetrized} divergence will be returned.}
}
\details{
  The Kullback-Leibler divergence (or \emph{relative entropy}) of two
  probability distributions \eqn{p} and \eqn{q} is defined as the
  integral
  \deqn{D_{\mathrm{KL}}(p\,||\,q) = \int_\Theta
    \log\Bigl(\frac{p(\theta)}{q(\theta)}\Bigr)\, p(\theta)\,
    \mathrm{d}\theta.}{D[KL](p || q) = Integral(log(p(theta)/q(theta))p(theta) d theta.}

  In the case of two normal distributions with mean and variance
  parameters given by (\eqn{\mu_1}{mu[1]}, \eqn{\Sigma_1}{Sigma[1]}) and
  (\eqn{\mu_2}{mu[2]}, \eqn{\Sigma_2}{Sigma[2]}), respectively, this
  results as 
  \deqn{D_{\mathrm{KL}}\bigl(p(\theta|\mu_1,\Sigma_1)\,||\,p(\theta|\mu_2,\Sigma_2)\bigr) = \frac{1}{2}\biggl(\mathrm{tr}(\Sigma_2^{-1} \Sigma_1) + (\mu_1-\mu_2)^\prime \Sigma_2^{-1} (\mu_1-\mu_2) - d + \log\Bigl(\frac{\det(\Sigma_2)}{\det(\Sigma_1)}\Bigr)\biggr)}{D[KL]( p(theta | mu[1],Sigma[1]) || p(theta | mu[2], Sigma[2]) ) = 0.5*(tr(Sigma[2]^-1 Sigma[1]) + (mu[1]-mu[2])' Sigma[2]^-1 (mu[1]-mu[2]) - d + log(det(Sigma[2]) / det(Sigma[1])))}
  where \eqn{d} is the dimension.

  The \emph{symmetrized} divergence simply results as
  \deqn{D_{\mathrm{s}}(p\,||\,q)=D_{\mathrm{KL}}(p\,||\,q)+D_{\mathrm{KL}}(q\,||\,p).}{D[s](p || q) = D[KL](p || q) + D[KL](q || p).}
}
\value{
  The divergence (\eqn{D_{\mathrm{KL}} \geq 0 }{D[KL] >= 0} or \eqn{D_{\mathrm{s}} \geq 0 }{D[s] >= 0}).
}
\references{
  S. Kullback. \emph{Information theory and statistics}.
  John Wiley and Sons, New York, 1959.

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
  \code{\link{bmr}}.
}
\examples{
kldiv(mu1=c(0,0), mu2=c(1,1), sigma1=diag(c(2,2)), sigma2=diag(c(3,3)))
}
\keyword{ distribution }
