\name{RhodesEtAlPrior}
\alias{RhodesEtAlPrior}
\alias{RhodesEtAlParameters}
\title{Heterogeneity priors for continuous outcomes (standardized mean
  differences) as proposed by Rhodes et al. (2015).}
\description{
  Use the prior specifications proposed in the paper by Rhodes et al.,
  based on an analysis of studies using standardized mean differences
  (SMD) that were published in the \emph{Cochrane Database of Systematic Reviews}.
}
\usage{
RhodesEtAlPrior(outcome=c(NA,"obstetric outcome",
  "resource use and hospital stay / process",
  "internal and external structure-related outcome",
  "general physical health and adverse event and pain and quality of life / functioning",
  paste("signs / symptoms reflecting continuation / end of condition and infection",
        "/ onset of new acute / chronic disease"),
  "mental health outcome", "biological marker", "various subjectively measured outcomes"),
  comparator1=c("pharmacological", "non-pharmacological", "placebo / control"),
  comparator2=c("pharmacological", "non-pharmacological", "placebo / control"),
  area=c("other","respiratory","cancer"))
}
\arguments{
  \item{outcome}{The type of outcome investigated (see below for a list
    of possible values). The default (\code{NA}) is the general
    (marginal) setting, without considering meta-analysis
    characteristics as covariates.}
  \item{comparator1}{One comparator's type.}
  \item{comparator2}{The other comparator's type.}
  \item{area}{The medical area.}
}
\details{ 
  Rhodes et al. conducted an analysis of studies listed in the
  \emph{Cochrane Database of Systematic Reviews} that were investigating
  standardized mean differences (SMD) as endpoints. As a result, they
  proposed empirically motivated log-Student-\eqn{t} prior distributions
  for the (squared!) heterogeneity parameter \eqn{\tau^2}, depending on
  the particular type of outcome investigated and the type of comparison
  in question. The underlying \eqn{t}-distribution's location and scale
  parameters here are internally stored in a 3-dimensional array (named
  \code{RhodesEtAlParameters}) and are most conveniently accessed using
  the \code{RhodesEtAlPrior()} function.
  
  The \code{outcome} argument specifies the type of outcome
  investigated. It may take one of the following values
  (partial matching is supported):
  \itemize{
    \item \code{NA}
    \item \code{"obstetric outcomes"}
    \item \code{"resource use and hospital stay / process"}
    \item \code{"internal and external structure-related outcome"}
    \item \code{"general physical health and adverse event and pain and quality of life / functioning"}
    \item \code{"signs / symptoms reflecting continuation / end of condition and infection / onset of new acute / chronic disease"}
    \item \code{"mental health outcome"}
    \item \code{"biological marker"}
    \item \code{"various subjectively measured outcomes"}.
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
    \item \code{"placebo / control"}.
  }
  Any combination is allowed for the \code{comparator1} and
  \code{comparator2} arguments, as long as not both arguments are set to
  \code{"placebo / control"}.
  The \code{area} argument specifies the medical context; possible
  values are:
  \itemize{
    \item \code{"respiratory"}
    \item \code{"cancer"}
    \item \code{"other"} (the default).
  }

  \bold{Note} that the location and scale parameters refer to the
  logarithmic (\emph{squared}) heterogeneity parameter \eqn{\tau^2},
  which is modelled using a Student-\eqn{t} distribution with 5 degrees
  of freedom.
  
}
\value{
  a list with elements
  \item{parameters}{the location and scale parameters (corresponding to the
    logarithmic, \emph{squared} heterogeneity parameter \eqn{\tau^2}).}
  \item{outcome.type}{the corresponding type of outcome.}
  \item{comparison.type}{the corresponding type of comparison.}
  \item{medical.area}{the medical context.}
  \item{dprior}{a \code{function(tau)} returning the prior
    density of \eqn{\tau}.}
  \item{pprior}{a \code{function(tau)} returning the prior cumulative
    distribution function (CDF) of \eqn{\tau}.}
  \item{qprior}{a \code{function(p)} returning the prior quantile
    function (inverse CDF) of \eqn{\tau}.}
}
\references{
  K.M. Rhodes, R.M. Turner, J.P.T. Higgins.
  \href{http://dx.doi.org/10.1016/j.jclinepi.2014.08.012}{Predictive
    distributions were developed for the extent of heterogeneity in
    meta-analyses of continuous outcome data}. \emph{Journal of Clinical
    Epidemiology}, \bold{68}(1):52-60, 2015.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{TurnerEtAlPrior}}.
}
\examples{
# determine prior distribution for a specific setting:
RP <- RhodesEtAlPrior("obstetric", "pharma", "placebo")
print(RP$parameters)
str(RP)
# a prior 95 percent interval for tau:
RP$qprior(c(0.025,0.975))

# the general (marginal) setting:
RP <- RhodesEtAlPrior()
print(RP$parameters)
str(RP)
# a prior 95 percent interval for tau:
RP$qprior(c(0.025,0.975))

\dontrun{
# load "metafor" package:
require("metafor")
# load data:
data("dat.normand1999")
# compute effect sizes (standardized mean differences):
es <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i,
             m2i=m2i, sd2i=sd2i, n2i=n2i,
             slab=source, data=dat.normand1999)

# derive appropriate prior:
RP <- RhodesEtAlPrior("resource use", "non-pharma", "non-pharma")
# show (central) prior 95 percent interval:
RP$qprior(c(0.025, 0.975))
# show prior 95 percent upper limit:
RP$qprior(0.95)

# perform meta analysis:
bma <- bayesmeta(es, tau.prior=RP$dprior)
# show results:
print(bma)
plot(bma, which=4, prior=TRUE)
}
}
\keyword{ distribution }
