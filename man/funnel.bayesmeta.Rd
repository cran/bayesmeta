\name{funnel.bayesmeta}
\alias{funnel.bayesmeta}
\title{
  Generate a funnel plot for a \code{\link{bayesmeta}} object.
}
\description{
  Generates a funnel plot, showing effect estimates (\eqn{y_i}{y[i]})
  vs. their standard errors (\eqn{\sigma_i}{sigma[i]}). 
}
\usage{
  \method{funnel}{bayesmeta}(x, main=deparse(substitute(x)), xlab=expression("effect "*y[i]),
       ylab=expression("standard error "*sigma[i]),
       zero=0.0, FE=FALSE, legend=FE, ...)
}
\arguments{
  \item{x}{
    a \code{\link{bayesmeta}} object.
  }
  \item{main}{
    main title for the plot.
  }
  \item{xlab}{
    x-axis title.
  }
  \item{ylab}{
    y-axis title.
  }
  \item{zero}{
    value at which a vertical \sQuote{reference} line should be drawn
    (default is 0). The line can be suppressed by setting this argument
    to \sQuote{\code{NA}}.
  }
  \item{FE}{
    a (\code{logical}) flag indicating whether the \dQuote{fixed effect}
    (FE) funnel (for \eqn{\tau=0}) is supposed to be shown along with
    the \dQuote{random effects} (RE) funnel.
  }
  \item{legend}{
    a (\code{logical}) flag indicating whether a legend identifying
    \sQuote{RE} and \sQuote{FE} funnels is supposed to be shown.
  }
  \item{...}{
    other arguments passed to the \code{plot()} function.
  }
}
\details{
  Generates a funnel plot of effect estimates (\eqn{y_i}{y[i]})
  on the x-axis vs. their associated standard errors
  (\eqn{\sigma_i}{sigma[i]}) on the y-axis (Note that the y-axis is
  pointing downwards). For many studies (large \eqn{k}) and in the
  absence of publication (selection) bias, the plot should resemble a
  (more or less) symmetric \dQuote{funnel} shape (Sterne \emph{et al},
  2005). Presence of publication bias, i.e., selection bias due to the
  fact that more dramatic effects may have higher chances of publication
  than less pronounced (or controversial) findings, may cause asymmetry
  in the plot; especially towards the bottom of the plot, studies may
  then populate a preferred corner. 

  Besides the \eqn{k} individual studies that are shown as circles, a
  vertical reference line is shown; its position is determined by
  the \sQuote{\code{zero}} argument.  The \dQuote{funnel} indicated in
  grey shows the estimated central 95\% prediction interval for
  \dQuote{new} effect estimates \eqn{y_i}{y[i]} conditional on a
  particular standard error \eqn{\sigma_i}{sigma[i]}, which results from
  convolving the prediction interval for the \emph{true} value
  \eqn{\theta_i}{theta[i]} with a normal distribution with variance
  \eqn{\sigma_i^2}{sigma[i]^2}. At \eqn{\sigma_i=0}{sigma[i]=0} (top of
  the funnel), this simply corresponds to the \dQuote{plain} prediction
  interval for \eqn{\theta_i}{theta[i]}. Convolutions are computed using
  the algorithm described in Roever and Friede (2017).

  By setting the \dQuote{\code{FE=TRUE}} argument, one may request a
  \dQuote{fixed effect} (FE) funnel along with the \dQuote{random
    effects} (RE) funnel that is shown by default. The FE funnel is
  analogous to the RE funnel, except that it is based on
  \emph{homogeneity} (\eqn{\tau=0}).
}
\note{Especially for few studies (small \eqn{k}), the conclusions from a
  forest plot are usually not very obvious (Sterne \emph{et al}, 2001;
  Terrin \emph{et al.}, 2005). Publication bias often requires rather
  large sample sizes to become apparent; funnel plots should hence
  always be interpreted with caution. 
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\references{J.A.C. Sterne, B.J. Becker and M. Egger.
  \href{https://doi.org/10.1002/0470870168.ch5}{The funnel plot}.
  In: H.R. Rothstein, A.J. Sutton and M. Borenstein, eds.
  \emph{Publication bias in meta-analysis - prevention,
        assessment and adjustments}.
  Wiley and Sons, 2005 (Chapter 5).

  N. Terrin, C.H. Schmid and J. Lau.
  \href{https://doi.org/10.1016/j.jclinepi.2005.01.006}{In an empirical
  evaluation of the funnel plot, researchers could not visually identify
  publication bias}. \emph{Journal of Clinical Epidemiology},
  \bold{58}(9):894-901, 2005.

  C. Roever, T. Friede.
  \href{https://doi.org/10.1080/10618600.2016.1276840}{Discrete
  approximation of a mixture distribution via restricted divergence}.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{26}(1):217-222, 2017.

} 
\seealso{
  \code{\link{bayesmeta}}, \code{\link[metafor]{funnel}}
}
\examples{
data("dat.egger2001", package="metafor")
es <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i,
             slab=study, data=dat.egger2001)
\dontrun{
bm <- bayesmeta(es)
print(bm)
forestplot(bm)
funnel(bm, xlab="logarithmic odds ratio", ylab="standard error",
       main="Egger (2001) example data")
}
}
\keyword{ hplot }
