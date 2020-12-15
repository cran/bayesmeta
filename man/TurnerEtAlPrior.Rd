\name{TurnerEtAlPrior}
\alias{TurnerEtAlPrior}
\alias{TurnerEtAlParameters}
\title{(Log-Normal) heterogeneity priors for binary outcomes
  as proposed by Turner et al. (2015).}
\description{
  Use the prior specifications proposed in the paper by Turner et al.,
  based on an analysis of studies using binary endpoints that were
  published in the \emph{Cochrane Database of Systematic Reviews}.
}
\usage{
TurnerEtAlPrior(outcome=c(NA, "all-cause mortality", "obstetric outcomes",
  "cause-specific mortality / major morbidity event / composite (mortality or morbidity)",
  "resource use / hospital stay / process", "surgical / device related success / failure",
  "withdrawals / drop-outs", "internal / structure-related outcomes",
  "general physical health indicators", "adverse events",
  "infection / onset of new disease",
  "signs / symptoms reflecting continuation / end of condition", "pain",
  "quality of life / functioning (dichotomized)", "mental health indicators",
  "biological markers (dichotomized)", "subjective outcomes (various)"),
  comparator1=c("pharmacological", "non-pharmacological", "placebo / control"),
  comparator2=c("pharmacological", "non-pharmacological", "placebo / control"))
}
\arguments{
  \item{outcome}{The type of outcome investigated (see below for a list
    of possible values).}
  \item{comparator1}{One comparator's type.}
  \item{comparator2}{The other comparator's type.}
}
\details{ 
  Turner et al. conducted an analysis of studies listed in the
  \emph{Cochrane Database of Systematic Reviews} that were investigating
  binary endpoints. As a result, they proposed empirically motivated
  log-normal prior distributions for the (squared!) heterogeneity
  parameter \eqn{\tau^2}, depending on the particular type of outcome
  investigated and the type of comparison in question. The log-normal
  parameters (\eqn{\mu} and \eqn{\sigma}) here are internally stored in
  a 3-dimensional array (named \code{TurnerEtAlParameters}) and are most
  conveniently accessed using the \code{TurnerEtAlPrior()} function.
  
  The \code{outcome} argument specifies the type of outcome
  investigated. It may take one of the following values
  (partial matching is supported):
  \itemize{
    \item \code{NA}
    \item \code{"all-cause mortality"}
    \item \code{"obstetric outcomes"}
    \item \code{"cause-specific mortality / major morbidity event / composite (mortality or morbidity)"}
    \item \code{"resource use / hospital stay / process"}
    \item \code{"surgical / device related success / failure"}
    \item \code{"withdrawals / drop-outs"}
    \item \code{"internal / structure-related outcomes"}
    \item \code{"general physical health indicators"}
    \item \code{"adverse events"}
    \item \code{"infection / onset of new disease"}
    \item \code{"signs / symptoms reflecting continuation / end of condition"}
    \item \code{"pain"}
    \item \code{"quality of life / functioning (dichotomized)"}
    \item \code{"mental health indicators"}
    \item \code{"biological markers (dichotomized)"}
    \item \code{"subjective outcomes (various)"}
  }
  Specifying \dQuote{\code{outcome=NA}} (the default) yields the
  \emph{marginal} setting, without considering meta-analysis
  characteristics as covariates.
  
  The \code{comparator1} and \code{comparator2} arguments together
  specify the type of comparison in question.  These may take one of the
  following values (partial matching is supported):
  \itemize{
    \item \code{"pharmacological"}
    \item \code{"non-pharmacological"}
    \item \code{"placebo / control"}
  }
  Any combination is allowed for the \code{comparator1} and
  \code{comparator2} arguments, as long as not both arguments are set to
  \code{"placebo / control"}.

  \bold{Note} that the log-normal prior parameters refer to the
  (\emph{squared}) heterogeneity parameter \eqn{\tau^2}. When you want
  to use the prior specifications for \eqn{\tau}, the square root,
  as the parameter (as is necessary when using the \code{bayesmeta()}
  function), you need to correct for the square root
  transformation. Taking the square root is equivalent to dividing by
  two on the log-scale, so the square root's distribution will still be
  log-normal, but with halved mean and standard deviation. The relevant
  transformations are already taken care of when using the resulting
  \code{$dprior()}, \code{$pprior()} and \code{$qprior()} functions; see
  also the example below.
}
\value{
  a list with elements
  \item{parameters}{the log-normal parameters (\eqn{\mu} and
    \eqn{\sigma}, corresponding to the \emph{squared} heterogeneity
    parameter \eqn{\tau^2} as well as \eqn{\tau}).}
  \item{outcome.type}{the corresponding type of outcome.}
  \item{comparison.type}{the corresponding type of comparison.}
  \item{dprior}{a \code{function(tau)} returning the prior
    density of \eqn{\tau}.}
  \item{pprior}{a \code{function(tau)} returning the prior cumulative
    distribution function (CDF) of \eqn{\tau}.}
  \item{qprior}{a \code{function(p)} returning the prior quantile
    function (inverse CDF) of \eqn{\tau}.}
}
\references{
  R.M. Turner, D. Jackson, Y. Wei, S.G. Thompson,
  J.P.T. Higgins. \href{https://doi.org/10.1002/sim.6381}{Predictive
    distributions for between-study heterogeneity and simple methods for
  their application in Bayesian meta-analysis}. \emph{Statistics in
  Medicine}, \bold{34}(6):984-998, 2015.

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
  \code{\link[stats]{dlnorm}}, \code{\link{RhodesEtAlPrior}}.
}
\examples{
# load example data:
data("CrinsEtAl2014")

# determine corresponding prior parameters:
TP <- TurnerEtAlPrior("surgical", "pharma", "placebo / control")
print(TP)
# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

\dontrun{
# compute effect sizes (log odds ratios) from count data
# (using "metafor" package's "escalc()" function):
crins.es <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)
print(crins.es)

# perform meta analysis:
crins.ma01 <- bayesmeta(crins.es, tau.prior=TP$dprior)
# for comparison perform analysis using weakly informative Cauchy prior:
crins.ma02 <- bayesmeta(crins.es, tau.prior=function(t){dhalfcauchy(t,scale=1)})

# show results:
print(crins.ma01)
print(crins.ma02)
# compare estimates; heterogeneity (tau):
rbind("Turner prior"=crins.ma01$summary[,"tau"], "Cauchy prior"=crins.ma02$summary[,"tau"])
# effect (mu):
rbind("Turner prior"=crins.ma01$summary[,"mu"], "Cauchy prior"=crins.ma02$summary[,"mu"])

# illustrate heterogeneity priors and posteriors:
par(mfcol=c(2,2))
  plot(crins.ma01, which=4, prior=TRUE, taulim=c(0,2),
       main="informative log-normal prior")
  plot(crins.ma02, which=4, prior=TRUE, taulim=c(0,2),
       main="weakly informative half-Cauchy prior")
  plot(crins.ma01, which=3, mulim=c(-3,0),
       main="informative log-normal prior")
  abline(v=0, lty=3)
  plot(crins.ma02, which=3, mulim=c(-3,0),
       main="weakly informative half-Cauchy prior")
  abline(v=0, lty=3)
par(mfrow=c(1,1))

# compare prior and posterior 95 percent upper limits for tau:
TP$qprior(0.95)
crins.ma01$qposterior(0.95)
qhalfcauchy(0.95)
crins.ma02$qposterior(0.95)
}
}
\keyword{ distribution }
