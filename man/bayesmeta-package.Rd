\name{bayesmeta-package}
\alias{bayesmeta-package}
\docType{package}
\title{
  Bayesian Random-Effects Meta-Analysis
}
\description{
  Description: A collection of functions allowing to derive the
  posterior distribution of the two parameters in a random-effects
  meta-analysis, and providing functionality to evaluate joint and
  marginal posterior probability distributions, predictive
  distributions, etc.
}
\details{
  \tabular{ll}{
    Package: \tab bayesmeta\cr
    Type:    \tab Package\cr
    Version: \tab 1.1\cr
    Date:    \tab 2015-12-18\cr
    License: \tab GPL (>=2)
  }
  The main functionality is provided by the \code{\link{bayesmeta}()}
  function. It takes the data (estimates and associated standard
  errors) and prior information (effect and heterogeneity priors), and
  returns an object containing functions that allow to derive posterior
  quantities like joint or marginal densities, quantiles, etc.
}
\author{
Christian Roever <christian.roever@med.uni-goettingen.de>
}
\references{
  Roever, C. et al. (2015). Bayesian random-effects meta-analysis. (in preparation).
}
\keyword{ models }
\keyword{ package }
\examples{
# example data by Snedecor and Cochran:
data(SnedecorCochran)

\dontrun{
# analysis using improper uniform prior
# (may take a few seconds to compute!):
bma <- bayesmeta(y=SnedecorCochran[,"mean"],
                 sigma=sqrt(SnedecorCochran[,"var"]),
                 label=SnedecorCochran[,"no"])

# show some summary statistics:
print(bma)

# show some plots:
plot(bma)
}
}
