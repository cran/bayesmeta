\name{summary.bmr}
\alias{summary.bmr}
\alias{print.summary.bmr}
\title{
  Summarizing a \code{\link{bmr}} object).
}
\description{
  Summarizes a \code{bmr} object, and (potentially) computes
  means and predictions.
}
\usage{
  \method{summary}{bmr}(object, X.mean, X.prediction, ...)
}
\arguments{
  \item{object}{
    a \code{\link{bmr}} object.
  }
  \item{X.mean}{
    a regressor matrix (\eqn{X}) for effect estimates that are to be
    derived. The matrix' row names define the labels passed on to the
    results.
  }
  \item{X.prediction}{
    an optional regressor matrix (\eqn{X}) for predictions that are to
    be derived. The matrix' row names define the labels passed on to the
    results.
  }
  \item{\ldots}{
    other arguments.
  }
}
\details{
  Prints details of the supplied \code{bmr} oject.

  Specification of the (optional) \dQuote{\code{X.mean}} or
  \dQuote{\code{X.prediction}} arguments allows to request computation
  of mean estimates or predictions corresponding to the supplied
  regressor matrices. Estimates (mode, median, mean, standard deviation,
  and 95 percent CI) may be retrieved from the returned object's
  \dQuote{\code{mean}} or \dQuote{\code{prediction}} elements (see
  example below).
}
\value{
  A list (of class \code{summary.bmr}) containing the following elements:
  \item{bmr}{the supplied \code{bmr} object.}
  \item{call}{an object of class \code{call} giving the
    function call that generated the \code{summary.bmr} object.}
  \item{X.mean, X.prediction}{the \sQuote{\code{X.mean}} and
    \sQuote{\code{X.prediction}} arguments.}
  \item{mean, prediction}{mean and predictions estimates (mode, median,
    mean, sd, and 95 percent credible intervals)}
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\references{
  C. Roever, T. Friede.
  Using the bayesmeta R package for Bayesian random-effects meta-regression.
  \emph{Computer Methods and Programs in Biomedicine},
  \bold{299}:107303, 2023.
  \doi{10.1016/j.cmpb.2022.107303}.
}
\examples{
\dontrun{
# perform a meta-analysis using binary ("indicator") covariables;
# load data:
data("CrinsEtAl2014")
# compute effect measures (log-OR):
crins.es <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)

# specify regressor matrix (binary indicator variables):
X <- cbind("basiliximab"=as.numeric(crins.es$IL2RA=="basiliximab"),
           "daclizumab" =as.numeric(crins.es$IL2RA=="daclizumab"))
print(X)

# perform meta-analysis:
bmr01 <- bmr(crins.es, X=X,
             tau.prior=function(t){dhalfnormal(t, scale=0.5)})

# show summary:
summary(bmr01)

# show summary with additional estimates and predictions:
summary(bmr01,
        X.mean = rbind("basiliximab" = c(1,0),
                       "daclizumab"  = c(0,1),
                       "difference"  = c(-1,1)),
        X.pred = rbind("basiliximab" = c(1,0),
                       "daclizumab"  = c(0,1)))

# compute mean estimates
smry <- summary(bmr01,
                X.mean = rbind("basiliximab" = c(1,0),
                               "daclizumab"  = c(0,1),
                               "difference"  = c(-1,1)))
# show mean estimates:
smry$mean
}
}
\keyword{ print }