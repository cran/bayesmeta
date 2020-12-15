\name{uisd}
\alias{uisd}
\alias{uisd.default}
\alias{uisd.escalc}
\concept{unit information standard deviation}
\title{
  Unit information standard deviation
}
\description{
  This function estimates the unit information standard deviation (UISD)
  from a given set of standard errors and associated sample sizes.
}
\usage{
  uisd(n, ...)
  \method{uisd}{default}(n, sigma, sigma2=sigma^2, labels=NULL, individual=FALSE, ...)
  \method{uisd}{escalc}(n, ...)
}
\arguments{
  \item{n}{
    vector of sample sizes \emph{or} an \code{\link[metafor]{escalc}} object.
}
  \item{sigma}{
    vector of standard errors associated with \code{n}.
}
  \item{sigma2}{
    vector of \emph{squared} standard errors (variances) associated with \code{n}.
}
  \item{labels}{
    (optional) a vector of labels corresponding to \code{n} and \code{sigma}.
}
  \item{individual}{
    a \code{logical} flag indicating whether individual (study-specific)
    UISDs are to be returned.
}
  \item{...}{
    other \code{uisd} arguments.
}
}
\details{
  The \emph{unit information standard deviation (UISD)} reflects the
  \dQuote{within-study} variability, which, depending on the effect
  measure considered, sometimes is a somewhat
  heuristic notion (Roever et al., 2020).
  For a single study, presuming that standard errors result as
  \deqn{\sigma_i=\frac{\sigma_\mathrm{u}}{\sqrt{n_i}},}{sigma[i] = sigma[u] / sqrt(n[i]),}
  where \eqn{\sigma_\mathrm{u}}{sigma[u]} is the within-study (population) standard
  deviation, the UISD simply results as
  \deqn{\sigma_\mathrm{u} = \sqrt{n_i \, \sigma_i^2}.}{sigma[u] = sqrt(n[i] * sigma[i]^2).}
  This is often appropriate when assuming an (approximately) normal likelihood.

  Assuming a constant \eqn{\sigma_\mathrm{u}}{sigma[u]} value across studies, this
  figure then may be estimated by
  \deqn{s_\mathrm{u} \;=\; \sqrt{\bar{n} \, \bar{s}^2_\mathrm{h}} \;=\; \sqrt{\frac{\sum_{i=1}^k n_i}{\sum_{i=1}^k \sigma_i^{-2}}},}{s[u] = sqrt(mean(n) * hmean(sigma^2)) = sqrt(sum(n)/sum(sigma^-2)),}
  where \eqn{\bar{n}}{mean(n)} is the average (arithmetic mean) of the
  studies' sample sizes, and \eqn{\bar{s}^2_\mathrm{h}}{hmean(sigma^2)} is the
  harmonic mean of the squared standard errors (variances).

  The estimator \eqn{s_\mathrm{u}}{s[u]} is motivated via meta-analysis
  using the normal-normal hierarchical model (NNHM). In the special case
  of homogeneity (zero heterogeneity, \eqn{\tau=0}{tau=0}), the overall
  mean estimate has standard error
  \deqn{\left(\sum_{i=1}^k\sigma_i^{-2}\right)^{-1/2}.}{sqrt(1/sum(sigma^(-2))).}
  Since this estimate corresponds to \emph{complete pooling}, the
  standard error may also be expressed via the UISD as
  \deqn{\frac{\sigma_\mathrm{u}}{\sqrt{\sum_{i=1}^k n_i}}.}{sigma[u] / sqrt(sum(n)).}
  Equating both above standard error expressions yields
  \eqn{s_\mathrm{u}}{s[u]} as an estimator
  of the UISD \eqn{\sigma_\mathrm{u}}{sigma[u]} (Roever \emph{et al}, 2020).
}
\value{
  Either a (single) estimate of the UISD, or, if \code{individual} was
  set to \sQuote{\code{TRUE}}, a (potentially named) vector of UISDs for
  each individual study.
}
\references{
  C. Roever, R. Bender, S. Dias, C.H. Schmid, H. Schmidli, S. Sturtz,
  S. Weber, T. Friede.
  On weakly informative prior distributions for the heterogeneity
  parameter in Bayesian random-effects meta-analysis.
  \emph{\href{https://arxiv.org/abs/2007.08352}{arXiv preprint 2007.08352}}
  (submitted for publication), 2020.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link[metafor]{escalc}}.
}
\examples{
# load data set:
data("CrinsEtAl2014")

# compute logarithmic odds ratios (log-ORs):
CrinsAR  <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)

# estimate the UISD:
uisd(n     = CrinsAR$exp.total + CrinsAR$cont.total,
     sigma = sqrt(CrinsAR$vi),
     label = CrinsAR$publication)

# for an "escalc" object, one may also apply the function directly:
uisd(CrinsAR)

# compute study-specific UISDs:
uisd(CrinsAR, individual=TRUE)
}
\keyword{ models }
