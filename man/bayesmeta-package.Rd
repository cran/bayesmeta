\name{bayesmeta-package}
\alias{bayesmeta-package}
\docType{package}
\title{
  Bayesian Random-Effects Meta-Analysis and Meta-Regression
}
\description{
  A collection of functions allowing to derive the posterior
  distribution of the model parameters in random-effects meta-analysis
  or meta-regression, and providing functionality to evaluate joint and
  marginal posterior probability distributions, predictive
  distributions, shrinkage effects, posterior predictive p-values, etc.
}
\details{
  \tabular{ll}{
    Package: \tab bayesmeta\cr
    Type:    \tab Package\cr
    Version: \tab 3.2\cr
    Date:    \tab 2022-09-29\cr
    License: \tab GPL (>=2)
  }
  The main functionality is provided by the \code{\link{bayesmeta}()}
  and \code{\link{bmr}()} function. It takes the data (estimates and
  associated standard errors) and prior information (effect and
  heterogeneity priors), and returns an object containing functions that
  allow to derive posterior quantities like joint or marginal densities,
  quantiles, etc. The \code{\link{bmr}()} function extends the approach
  to meta-regression by allowing to specify covariables (moderators) in
  addition.
}
\author{
Christian Roever <christian.roever@med.uni-goettingen.de>
}
\seealso{
  \code{\link{forestplot.bayesmeta}}, \code{\link{plot.bayesmeta}},
  \code{\link{bmr}}.
}
\references{
  C. Roever.
  Bayesian random-effects meta-analysis using the bayesmeta R package.
  \emph{Journal of Statistical Software}, \bold{93}(6):1-51, 2020.
  \doi{10.18637/jss.v093.i06}.
  
  C. Roever, T. Friede.
  Using the bayesmeta R package for Bayesian random-effects meta-regression.
  \href{https://arxiv.org/abs/2209.06004}{arXiv preprint 2209.06004}, 2022.
}
\keyword{ models }
\keyword{ package }
\examples{
# example data by Snedecor and Cochran:
data("SnedecorCochran")

\dontrun{
# analysis using improper uniform prior
# (may take a few seconds to compute!):
bma <- bayesmeta(y=SnedecorCochran[,"mean"],
                 sigma=sqrt(SnedecorCochran[,"var"]),
                 label=SnedecorCochran[,"no"])

# show some summary statistics:
bma

# show a bit more details:
summary(bma)

# show a forest plot:
forestplot(bma)

# show some more plots:
plot(bma)
}
}
