\name{forestplot.bayesmeta}
\alias{forestplot.bayesmeta}
\title{
  Generate a forest plot for a \code{\link{bayesmeta}} object
  (based on the \code{forestplot} package's plotting functions).
}
\description{
  Generates a forest plot, showing individual estimates along with their
  95 percent confidence intervals, shrinkage intervals, resulting effect
  estimate and prediction interval.
}
\usage{
  \method{forestplot}{bayesmeta}(x, labeltext, exponentiate=FALSE,
           prediction=TRUE, shrinkage=TRUE, heterogeneity=TRUE, 
           digits=2, plot=TRUE,
           fn.ci_norm, fn.ci_sum, col, legend, boxsize, ...)
}
\arguments{
  \item{x}{
    a \code{\link{bayesmeta}} object.
  }
  \item{labeltext}{an (alternative) \dQuote{\code{labeltext}} argument
    which is then handed on to the \code{\link[forestplot]{forestplot}()}
    function (see the help there). You can use this to change contents
    or add columns to the displayed table; see the example below.
  }
  \item{exponentiate}{
    a logical flag indicating whether to exponentiate numbers (effect
    sizes) in table and plot.
  }
  \item{prediction}{
    a logical flag indicating whether to show the prediction interval
    below the mean estimate.
  }
  \item{shrinkage}{
    a logical flag indicating whether to show shrinkage intervals along
    with the quoted estimates.
  }
  \item{heterogeneity}{
    a logical flag indicating whether to quote the heterogeneity estimate 
    and CI (at the bottom left).
  }
  \item{digits}{
    The number of significant digits to be shown.
    This is interpreted relative to the standard errors of all estimates.
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
  resulting estimates (effect estimate and prediction interval,
  as well as shrinkage estimates and intervals).
}
\note{This function is based on the \pkg{forestplot} package's
      \dQuote{\code{\link[forestplot]{forestplot}()}} function.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\references{
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
  \code{\link{forest.bayesmeta}},
  \code{\link{plot.bayesmeta}}.
}
\examples{
# load data:
data("CrinsEtAl2014")

\dontrun{
# compute effect sizes (log odds ratios) from count data
# (using "metafor" package's "escalc()" function):
require("metafor")
crins.es <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)
print(crins.es)

# perform meta analysis:
crins.ma <- bayesmeta(crins.es, tau.prior=function(t){dhalfcauchy(t,scale=1)})

########################
# generate forest plots
require("forestplot")

# default options:
forestplot(crins.ma)

# exponentiate values (shown in table and plot), show vertical line at OR=1:
forestplot(crins.ma, expo=TRUE, zero=1)

# logarithmic x-axis:
forestplot(crins.ma, expo=TRUE, xlog=TRUE)

# omit prediction interval:
forestplot(crins.ma, predict=FALSE)

# omit shrinkage intervals:
forestplot(crins.ma, shrink=FALSE)

# show more decimal places:
forestplot(crins.ma, digits=3)

# change table values:
# (here: add columns for event counts)
fp <- forestplot(crins.ma, expo=TRUE, plot=FALSE)
labtext <- fp$labeltext
labtext <- cbind(labtext[,1],
                 c("treatment",
                   paste0(CrinsEtAl2014[,"exp.AR.events"], "/", CrinsEtAl2014[,"exp.total"]),
                   "",""),
                 c("control",
                   paste0(CrinsEtAl2014[,"cont.AR.events"], "/", CrinsEtAl2014[,"cont.total"]),
                   "",""),
                 labtext[,2:3])
labtext[1,4] <- "OR"
print(fp$labeltext) # before
print(labtext)      # after
forestplot(crins.ma, labeltext=labtext, expo=TRUE, xlog=TRUE)

# see also the "forestplot" help for more arguments that you may change,
# e.g. the "clip", "xticks", "xlab" and "title" arguments,
# or the "txt_gp" argument for label sizes etc.:
forestplot(crins.ma, clip=c(-4,1), xticks=(-3):0,
           xlab="log-OR", title="pediatric transplantation example",
           txt_gp = fpTxtGp(ticks = gpar(cex=1), xlab = gpar(cex=1)))
}
}
\keyword{ hplot }
