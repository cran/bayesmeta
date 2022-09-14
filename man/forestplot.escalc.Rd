\name{forestplot.escalc}
\alias{forestplot.escalc}
\title{
  Generate a forest plot for an \code{\link[metafor]{escalc}} object
  (based on the \code{forestplot} package's plotting functions).
}
\description{
  Generates a forest plot, showing estimates along with their
  95 percent confidence intervals.
}
\usage{
  \method{forestplot}{escalc}(x, labeltext, exponentiate=FALSE,
           digits=2, plot=TRUE,
           fn.ci_norm, fn.ci_sum, col, legend, boxsize, ...)
}
\arguments{
  \item{x}{
    an \code{\link[metafor]{escalc}} object.
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
  Generates a forest plot illustrating the data returned by the
  \dQuote{\code{\link[metafor]{escalc}()}} function (showing the
  estimates potentially to be meta-analyzed, but without a combined
  summary estimate).
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

  R.D. Riley, J.P. Higgins and J.J. Deeks.
  Interpretation of random effects meta-analyses.
  \emph{BMJ}, \bold{342}:d549, 2011.
  \doi{10.1136/bmj.d549}.
} 
\seealso{
  \code{\link[metafor]{escalc}},
  \code{\link[forestplot]{forestplot}},
  \code{\link{forestplot.bayesmeta}},
  \code{\link{forestplot.bmr}}.
}
\examples{
# load data:
data("CrinsEtAl2014")

# compute effect sizes (log odds ratios) from count data
# (using "metafor" package's "escalc()" function):
crins.es <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)
print(crins.es)

########################
# generate forest plots;
# with default settings:
forestplot(crins.es)

# exponentiate values (shown in table and plot), show vertical line at OR=1:
forestplot(crins.es, expo=TRUE, zero=1)

# logarithmic x-axis:
forestplot(crins.es, expo=TRUE, xlog=TRUE)

# show more decimal places:
forestplot(crins.es, digits=3)

# change table values:
# (here: add columns for event counts)
fp <- forestplot(crins.es, expo=TRUE, plot=FALSE)
labtext <- fp$labeltext
labtext <- cbind(labtext[,1],
                 c("treatment",
                   paste0(CrinsEtAl2014[,"exp.AR.events"], "/", CrinsEtAl2014[,"exp.total"])),
                 c("control",
                   paste0(CrinsEtAl2014[,"cont.AR.events"], "/", CrinsEtAl2014[,"cont.total"])),
                 labtext[,2:3])
labtext[1,4] <- "OR"
print(fp$labeltext) # before
print(labtext)      # after
forestplot(crins.es, labeltext=labtext, expo=TRUE, xlog=TRUE)

# see also the "forestplot" help for more arguments that you may change,
# e.g. the "clip", "xticks", "xlab" and "title" arguments,
# or the "txt_gp" argument for label sizes etc.:
forestplot(crins.es, clip=c(-4,1), xticks=(-3):0,
           xlab="log-OR", title="pediatric transplantation example",
           txt_gp = fpTxtGp(ticks = gpar(cex=1), xlab = gpar(cex=1)))


###########################################################
# In case effects and standard errors are computed already
# (and normally one wouldn't need to call "escalc()")
# you can still use "escalc()" to assemble the plot, e.g.:

data("HinksEtAl2010")
print(HinksEtAl2010)

hinks.es <- escalc(yi=log.or, sei=log.or.se,
                   slab=study, measure="OR",
                   data=HinksEtAl2010)

forestplot(hinks.es)
}
\keyword{ hplot }
