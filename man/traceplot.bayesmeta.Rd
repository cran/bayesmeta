\name{traceplot}
\alias{traceplot}
\alias{traceplot.default}
\alias{traceplot.bayesmeta}
\alias{traceplot.bmr}
\title{
  Illustrate conditional means of study-specific estimates as well as
  overall mean (or other linear combinations) as a function of
  heterogeneity.
}
\description{
  Generates a trace plot of study-specific (shrinkage) estimates as a
  function of the heterogeneity (\eqn{\tau}), along with conditional
  estimates of the overall mean or other linear combinations of
  regression parameters.
  The heterogeneity's posterior distribution is also shown at the bottom.
}
\usage{
  traceplot(x, ...)
  \method{traceplot}{bayesmeta}(x, mulim, taulim, ci=FALSE,
          ylab="effect", prior=FALSE, infinity=FALSE,
          rightmargin=8, col=rainbow(x$k), labcol=col,
          meanlabel="overall mean", meancol="black",
          meanlabcol=meancol, ...)
  \method{traceplot}{bmr}(x, mulim, taulim, ci=FALSE,
          ylab="effect", prior=FALSE, infinity=FALSE,
          rightmargin=8, col=rainbow(x$k), labcol=col,
          X, Xlabels, Xcols="black", Xlabcols=Xcols, ...)
}
\arguments{
  \item{x}{
    a \code{\link{bayesmeta}} or \code{\link{bmr}}object.
  }
  \item{mulim, taulim}{(optional) ranges for the effect (mu) and
    heterogeneity (tau) axes. If only one value is given for
    \code{taulim}, then this is taken as the upper limit, and the lower
    limit is zero.
  }
  \item{ci}{
    a logical flag indicating whether to also show (conditional)
    confidence intervals.
  }
  \item{ylab}{
    a label for the effect (mu) axis.
  }
  \item{prior}{
    a logical flag indicating whether to show the (heterogeneity) prior
    density along with its posterior.
  }
  \item{infinity}{
    a logical flag indicating whether add an \dQuote{infinity} tickmark
    to the heterogeneity (tau) axis and show the corresponding limiting
    values.
  }
  \item{rightmargin}{
    an additional margin to be added to the right side of the plot, in
    order to accomodate the estimates' labels. In case study labels
    still extend beyond the figure margin, try increasing this number.
  }
  \item{col}{
    colors to be used for plotting the (\eqn{k}) estimates.
  }
  \item{labcol}{
    colors to be used for labelling the (\eqn{k}) estimates.
  }
  \item{meanlabel}{
    a label for the overall mean estimate
    (\code{traceplot.bayesmeta()}).
  }
  \item{meancol}{
    colour specification for the overall mean estimate
    (\code{traceplot.bayesmeta()}).
  }
  \item{meanlabcol}{
    colour specification for the overall mean label
    (\code{traceplot.bayesmeta()}).
  }
  \item{X}{
    matrix (or vector) of coefficients defining linear combinations of
    regression parameters to be shown (\code{traceplot.bmr()}).
  }
  \item{Xlabels}{
    labels for the linear combinations (\code{traceplot.bmr()}).
  }
  \item{Xcols}{
    colour specification for the linear combinations
    (\code{traceplot.bmr()}).
  }
  \item{Xlabcols}{
    colour specification for the linear combination labels
    (\code{traceplot.bmr()}).
  }
  \item{...}{
    other arguments passed on to the
    \code{\link[graphics]{plot}()} function.
  }
}
\details{
  For a given heterogeneity (\eqn{\tau}) value, the \emph{conditional}
  posterior distributions of the overall effect (\eqn{\mu}) as well as
  the study-specific parameters (\eqn{\theta_i}{theta[i]}) are again
  normal. The conditional normal moments (mean and variance) then vary
  as functions of the heterogeneity; for large heterogeneity, the
  shrinkage estimates approach the original data (\eqn{y_i}{y[i]}),
  while the overall mean approaches an un-weighted overall average. For
  small heterogeneity, both overall mean as well as study-specific
  estimates are increasingly \emph{shrunk} towards the
  inverse-variance-weighted \sQuote{common-effect} estimate (Roever,
  2020).

  This trace plot illustrates the conditional (overall and
  study-specific) estimates along with the heterogeneity's posterior
  distribution (density) in a layout similar to that utilized
  by Rubin (1981).
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\references{
  C. Roever.
  Bayesian random-effects meta-analysis using the bayesmeta R package.
  \emph{Journal of Statistical Software}, \bold{93}(6):1-51, 2020.
  \doi{10.18637/jss.v093.i06}.
  
  C. Roever, T. Friede.
  Using the bayesmeta R package for Bayesian random-effects meta-regression.
  \emph{Computer Methods and Programs in Biomedicine},
  \bold{299}:107303, 2023.
  \doi{10.1016/j.cmpb.2022.107303}.

  D.B. Rubin.
  Estimation in parallel randomized experiments.
  \emph{Journal of Educational Statistics}, \bold{6}(4):377-401, 1981.
  \doi{10.3102/10769986006004377}.
  
  DuMouchel, W. H. (1994).
  Hierarchical Bayes linear models for meta-analysis.
  Technical Report 27, National Institute of Statistical Sciences (NISS);
  Research Triangle Park, NC, USA.
  \url{https://www.niss.org/research/technical-reports/hierarchical-bayes-linear-models-meta-analysis-1994}
} 
\seealso{
  \code{\link{bayesmeta}}, \code{\link{bmr}}.
}
\examples{
\dontrun{
########################
# SAT coaching example;
# load example data:
data("Rubin1981")

# perform meta-analysis:
bma01 <- bayesmeta(y=Rubin1981[,"effect"], sigma=Rubin1981[,"stderr"],
                  labels=Rubin1981[,"school"], tau.prior="uniform")

# show meta-analysis results:
forestplot(bma01)

# show trace plot:
traceplot(bma01)


##################################
# COPD (meta-regression) example;
# load example data,
# compute effect sizes (log-ORs):
data("KarnerEtAl2014")
karner.exa <- escalc(measure="OR",
                     ai=tiotropium.exa, n1i=tiotropium.total,
                     ci=placebo.exa,    n2i=placebo.total,
                     slab=study, data=KarnerEtAl2014)


#################################
# perform "plain" meta-analysis:
bma02 <- bayesmeta(karner.exa,
                   tau.prior=function(t){dhalfnormal(t, scale=0.5)})

traceplot(bma02, ylab="log-OR",
          prior=TRUE, infi=TRUE, taulim=0.53)

forestplot(bma02)


############################
# perform meta-regressions:

# 1st regression;
# specify regressor matrix
# (indicator variables, "short" vs. "long" study duration):
X1 <- cbind("short" = as.numeric(karner.exa$duration == "up to 1 year"),
            "long"  = as.numeric(karner.exa$duration == "1 year or longer"))

# perform meta-regression
# (two group means, common heterogeneity):
bmr01 <- bmr(karner.exa, X=X1,
             tau.prior=function(t){dhalfnormal(t, scale=0.5)})

# show trace plot:
traceplot(bmr01, ylab="log-OR", prior=TRUE,
          taulim=0.53, mulim=c(-1, 0.2),
          X=rbind("short" = c(1,0),
                  "long"  = c(0,1)))

# 2nd regression;
# specify regressor matrix
# (baseline FEV1, an indicator of disease severity):
X2 <- cbind("intercept" = 1,
            "fev1"      = karner.exa$baseline.fev1)

# perform meta-regression
# (linear effect of FEV1 on log-OR):
bmr02 <- bmr(karner.exa, X=X2,
             tau.prior=function(t){dhalfnormal(t, scale=0.5)})

traceplot(bmr02, ylab="log-OR", prior=TRUE,
          taulim=0.53, mulim=c(-1.0, 0.2),
          X=rbind("FEV1 = 1.0"=c(1,1.0),
                  "FEV1 = 1.5"=c(1,1.5),
                  "FEV1 = 2.0"=c(1,2.0)))
}
}
\keyword{ hplot }
