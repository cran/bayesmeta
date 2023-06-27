\name{forestplot.bmr}
\alias{forestplot.bmr}
\title{
  Generate a forest plot for a \code{\link{bmr}} object
  (based on the \code{forestplot} package's plotting functions).
}
\description{
  Generates a forest plot, showing individual estimates along with their
  95 percent confidence intervals, shrinkage intervals, resulting effect
  estimates and prediction intervals.
}
\usage{
  \method{forestplot}{bmr}(x, X.mean, X.prediction,
           labeltext, exponentiate=FALSE,
           shrinkage=TRUE, heterogeneity=TRUE, 
           digits=2, decplaces.X, plot=TRUE,
           fn.ci_norm, fn.ci_sum, col, legend=NULL, boxsize, ...)
}
\arguments{
  \item{x}{
    a \code{\link{bmr}} object.
  }
  \item{X.mean}{
    a regressor matrix (\eqn{X}) for effect estimates that are to be
    shown in the plot. By default, a diagonal matrix; set to
    \code{NULL} in order to suppress showing summary estimates. The
    matrix' row names define the labels shown in the plot.
  }
  \item{X.prediction}{
    an optional regressor matrix (\eqn{X}) for predictions that are to
    be shown in the plot. The matrix' row names define the labels shown
    in the plot.
  }
  \item{labeltext}{an (alternative) \dQuote{\code{labeltext}} argument
    which is then handed on to the \code{\link[forestplot]{forestplot}()}
    function (see the help there). You can use this to change contents
    or add columns to the displayed table; see also the example below.
  }
  \item{exponentiate}{
    a logical flag indicating whether to exponentiate numbers (effect
    sizes) in table and plot.
  }
  \item{shrinkage}{
    a logical flag indicating whether to show shrinkage intervals along
    with the quoted estimates.
  }
  \item{heterogeneity}{
    a logical flag indicating whether to quote the heterogeneity estimate 
    and CI (at the bottom left of the plot).
  }
  \item{digits}{
    the number of significant digits to be shown.
    This is interpreted relative to the standard errors of all estimates.
  }
  \item{decplaces.X}{
    The number of decimal places to be shown for the regressors.
  }
  \item{plot}{
    a logical flag indicating whether to actually generate a plot.
  }
  \item{fn.ci_norm, fn.ci_sum, col, legend, boxsize, \ldots}{
    other arguments passed on to the
    \pkg{forestplot} package's \code{\link[forestplot]{forestplot}}
    function (see also the help there).
  }
}
\details{
  Generates a forest plot illustrating the underlying data and
  resulting estimates (effect estimates and/or prediction intervals,
  as well as shrinkage estimates and intervals).
  For effect estimates and prediction intervals, regressor matrices
  (\eqn{X}) need to be supplied via the \sQuote{\code{X.mean}} or
  \sQuote{\code{X.prediction}} arguments. Effect estimates are shown as
  diamonds, predictions are shown as horizontal bars.
}
\note{This function is based on the \pkg{forestplot} package's
      \dQuote{\code{\link[forestplot]{forestplot}()}} function.
}
\value{
  A list containing the following elements:
  \item{data}{a \code{matrix} of estimates and CIs.}
  \item{X.mean, X.prediction}{the \sQuote{\code{X.mean}} and
    \sQuote{\code{X.prediction}} arguments.}
  \item{shrinkage}{a \code{matrix} of shrinkage estimates and CIs.}
  \item{labeltext}{a \code{matrix} of table entries.}
  \item{forestplot}{result of the call to the
    \sQuote{\code{forestplot()}} function.}
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

  C. Roever.
  Bayesian random-effects meta-analysis using the bayesmeta R package.
  \emph{Journal of Statistical Software}, \bold{93}(6):1-51, 2020.
  \doi{10.18637/jss.v093.i06}.

  C. Lewis and M. Clarke.
  Forest plots: trying to see the wood and the trees.
  \emph{BMJ}, \bold{322}:1479, 2001.
  \doi{10.1136/bmj.322.7300.1479}.

  C. Guddat, U. Grouven, R. Bender and G. Skipka.
  A note on the graphical presentation of prediction intervals in
  random-effects meta-analyses.
  \emph{Systematic Reviews}, \bold{1}(34), 2012.
  \doi{10.1186/2046-4053-1-34}.
  
  R.D. Riley, J.P. Higgins and J.J. Deeks.
  Interpretation of random effects meta-analyses.
  \emph{BMJ}, \bold{342}:d549, 2011.
  \doi{10.1136/bmj.d549}.
} 
\seealso{
  \code{\link{bayesmeta}},
  \code{\link[forestplot]{forestplot}},
  \code{\link{forestplot.bayesmeta}},
  \code{\link{forestplot.escalc}}.
}
\examples{
\dontrun{
#################################################################
# perform a meta-analysis using binary ("indicator") covariables:

# load data:
data("CrinsEtAl2014")
# compute effect measures (log-OR):
crins.es <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)
# show data:
crins.es[,c("publication", "IL2RA", "exp.AR.events", "exp.total",
            "cont.AR.events", "cont.total", "yi", "vi")]
# specify regressor matrix (binary indicator variables):
X <- cbind("basiliximab"=as.numeric(crins.es$IL2RA=="basiliximab"),
           "daclizumab" =as.numeric(crins.es$IL2RA=="daclizumab"))
print(X)
# perform meta-analysis:
bmr01 <- bmr(crins.es, X=X,
             tau.prior=function(t){dhalfnormal(t, scale=0.5)})

# show forest plot:
forestplot(bmr01)

# show forest plot including contrast
# (difference between the two groups):
forestplot(bmr01,
           X.mean=rbind("basiliximab"      = c(1, 0),
                        "daclizumab"       = c(0, 1),
                        "group difference" = c(-1, 1)))

##############################################
# perform the meta-analysis using a different 
# ("intercept / slope") regressor setup:
X <- cbind("intercept"=1,
           "offset.dac"=as.numeric(crins.es$IL2RA=="daclizumab"))
print(X)
# perform meta-analysis:
bmr02 <- bmr(crins.es, X=X,
             tau.prior=function(t){dhalfnormal(t, scale=0.5)})

# show default forest plot:
forestplot(bmr02)

# show forest plot including both group means and their difference:
forestplot(bmr02,
           X.mean=rbind("basiliximab"      = c(1, 0),
                        "daclizumab"       = c(1, 1),
                        "group difference" = c(0, 1)))

###############################################################
# a meta analysis using a continuous regressor
# and including prediction:
help("NicholasEtAl2019")
# load data:
data("NicholasEtAl2019")
# compute effect sizes (logarithic odds) from count data:
es <- escalc(measure="PLO",
             xi=patients*(prog.percent/100), ni=patients,
             slab=study, data=NicholasEtAl2019)
# set up regressor matrix:
X <- cbind("intercept2000" = 1, "year" = (es$year-2000))
# perform analysis:
bmr03 <- bmr(es, X=X)
# show forest plot including some mean estimates for the
# years from 1990 to 2018, and a prediction for 2019:
forestplot(bmr03,
           X.mean=rbind("intercept (2000)"  = c(1,   0),
                        "annual change"     = c(0,   1),
                        "change per decade" = c(0,  10),
                        "mean 1990"         = c(1, -10),
                        "mean 2000"         = c(1,   0),
                        "mean 2010"         = c(1,  10),
                        "mean 2018"         = c(1,  18)),
           X.predict=rbind("prediction 2019"   = c(1,  19)),
           xlab="log-odds",
           txt_gp = fpTxtGp(ticks = gpar(cex=1), xlab = gpar(cex=1)))

# the shown summaries and predictions may also be computed "manually";
# mean effect (year 2018), posterior median and 95 percent CI:
bmr03$qpredict(p=0.5, x=c(1, 18))
bmr03$pred.interval(level=0.95, x=c(1, 18))

# prediction (year 2019), posterior median and 95 percent CI:
bmr03$qpredict(p=0.5, x=c(1, 19), mean=FALSE)
bmr03$pred.interval(level=0.95, x=c(1, 19), mean=FALSE)

# means and predictions may also be derived
# using the "summary()" function:
summary(bmr03,
        X.mean=rbind("intercept (2000)"  = c(1,   0),
                     "annual change"     = c(0,   1),
                     "change per decade" = c(0,  10),
                     "mean 1990"         = c(1, -10),
                     "mean 2000"         = c(1,   0),
                     "mean 2010"         = c(1,  10),
                     "mean 2018"         = c(1,  18)),
        X.predict=rbind("prediction 2019"   = c(1,  19)))


##########################################################
# the tabular part of the forest plot may also be changed;
# draw a default plot:
forestplot(bmr03)
# don't plot, only extract the tabular bits:
fp <- forestplot(bmr03, plot=FALSE)
labtxt <- fp$labeltext
head(labtxt)

# drop two columns:
labtxt <- labtxt[,-c(2,3)]
# add two new columns:
labtxt <- cbind(labtxt[,1],
                c("year", es$year, "",""),
                c("events / total",
                  paste(round(es$patients*(es$prog.percent/100)),
                  "/", es$patients), "",""),
                labtxt[,2:3])
head(labtxt)
# draw new forest plot:
forestplot(bmr03, labeltext=labtxt, xlab="log-odds")
}
}
\keyword{ hplot }