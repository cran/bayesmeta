\name{ess}
\alias{ess}
\alias{ess.bayesmeta}
\concept{Effective sample size (ESS)}
\title{
  Effective sample size (ESS)
}
\description{
  This function computes the effective sample size (ESS) of a posterior
  predictive distribution.
}
\usage{
  ess(object, ...)
  \method{ess}{bayesmeta}(object, uisd, method=c("elir", "vr", "pr", "mtm.pt"), ...)
}
\arguments{
  \item{object}{a \code{bayesmeta} object.}
  \item{uisd}{the \emph{unit infomation standard deviation}
    (a single numerical value, \emph{or} a \code{function} of the
    parameter (\eqn{\mu}{mu})).}
  \item{method}{a character string specifying the method to be used for
    ESS computation. By default, the expected local-information-ratio
    ESS (\eqn{ESS_{ELIR}}{ESS_ELIR}) is returned.}
  \item{...}{additional arguments}
}
\details{
  The information conveyed by a prior distribution may often be
  quantified in terms of an \emph{effective sample size
    (ESS)}. Meta-analyses are commonly utilized to summarize
  \dQuote{historical} information in order to inform a future study,
  leading to a \emph{meta-analytic-predictive (MAP) prior} (Schmidli et
  al., 2014). In the context of the normal-normal hierarchical model
  (NNHM), the MAP prior results as the (posterior) predictive
  distribution for a \dQuote{new} study mean
  \eqn{\theta_{k+1}}{theta[k+1]}. This function computes the ESS for the
  posterior predictive distribution based on a \code{bayesmeta} object.

  Within the NNHM, the notion of an effective sample size requires the
  specification of a \emph{unit information standard deviation (UISD)}
  (Roever et al., 2020); see also the \sQuote{\code{\link{uisd}()}}
  function's help page. The UISD \eqn{\sigma_\mathrm{u}}{sigma[u]} here
  determines the \emph{Fisher information for one information unit},
  effectively assuming that a study's sample size \eqn{n_i}{n[i]} and
  standard error \eqn{\sigma_i}{sigma[i]} are related simply as
  \deqn{\sigma_i=\frac{\sigma_\mathrm{u}}{\sqrt{n_i}},}{sigma[i] = sigma[u] / sqrt(n[i]),}
  i.e., the squared standard error is inversely proportional to the
  sample size. For the (possibly hypothetical) case of a sample size of
  \eqn{n_i=1}{n[i]=1}, the standard error then is equal to the UISD
  \eqn{\sigma_\mathrm{u}}{sigma[u]}.

  Specifying the UISD as a constant is often an approximation,
  sometimes it is also possible to specify the UISD as a function of the
  parameter (\eqn{\mu}{mu}). For example, in case the outcome in
  the meta-analyses are log-odds, then the UISD varies with the (log-)
  odds and is given by \eqn{2\,\mathrm{cosh}(\mu/2)}{2*cosh(mu/2)}
  (see also the example below).

  The ESS may be computed or approximated in several ways.
  Possible choices here are:
  \itemize{
    \item \code{"elir"}: the \emph{expected local-information-ratio (ELIR)} method (the default),
    \item \code{"vr"}: the \emph{variance ratio (VR)} method,
    \item \code{"pr"}: the \emph{precision ratio (PR)} method,
    \item \code{"mtm.pt"}: the \emph{Morita-Thall-Mueller / Pennello-Thompson (MTM.PM)} method.
  }
  For more details on these see also Neuenschwander et al. (2020).
  
}
\value{
  The effective sample size (ESS).
}
\references{
  B. Neuenschwander, S. Weber, H. Schmidli, A. O'Hagan.
  Predictively consistent prior effective sample sizes.
  \emph{Biometrics}, \bold{76}(2):578-587, 2020.
  \doi{10.1111/biom.13252}.
  
  H. Schmidli, S. Gsteiger, S. Roychoudhury, A. O'Hagan,
  D. Spiegelhalter, B. Neuenschwander.
  Robust meta-analytic-predictive priors in clinical trials with
  historical control information.
  \emph{Biometrics}, \bold{70}(4):1023-1032, 2014.
  \doi{10.1111/biom.12242}.

  C. Roever, R. Bender, S. Dias, C.H. Schmid, H. Schmidli, S. Sturtz,
  S. Weber, T. Friede.
  On weakly informative prior distributions for the heterogeneity
  parameter in Bayesian random-effects meta-analysis.
  \emph{Research Synthesis Methods}, \bold{12}(4):448-474, 2021.
  \doi{10.1002/jrsm.1475}.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{bayesmeta}}, \code{\link{uisd}}.
}
\examples{
# load data set:
data("BaetenEtAl2013")
print(BaetenEtAl2013)

\dontrun{
# compute effect sizes (logarithmic odds) from the count data:
as <- escalc(xi=events, ni=total, slab=study,
             measure="PLO", data=BaetenEtAl2013)

# estimate the unit information standard deviation (UISD):
uisd(as, individual=TRUE)
uisd(as)  # = 2.35

# perform meta-analysis
# (using uniform priors for effect and heterogeneity):
bm <- bayesmeta(as)

# show forest plot:
forestplot(bm, zero=NA, xlab="log-odds")

# compute ESS_ELIR (based on fixed UISD):
ess(bm, uisd=2.35)  # = 45.7 patients

# compute ESS_ELIR based on UISD as a function of the log-odds:
uisdLogOdds <- function(logodds)
{
  return(2 * cosh(logodds / 2))
}

# Note: in the present example, probabilities are
# at approximately 0.25, corresponding to odds of 1:3.
uisdLogOdds(log(1/3))
# The UISD value of 2.31 roughly matches the above empirical figure.

ess(bm, uisd=uisdLogOdds)  # = 43.4 patients
}
}
\keyword{ models }
