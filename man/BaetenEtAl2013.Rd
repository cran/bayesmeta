\name{BaetenEtAl2013}
\docType{data}
\alias{BaetenEtAl2013}
\title{Ankylosing spondylitis example data}
\description{Numbers of cases (patients) and events (responders)
  in the placebo control groups of eight studies.}
\usage{data("BaetenEtAl2013")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{study}  \tab \code{character} \tab study label \cr
    \bold{year}   \tab \code{numeric}   \tab publication year \cr
    \bold{events} \tab \code{numeric}   \tab number of responders \cr
    \bold{total}  \tab \code{numeric}   \tab total number of patients
  }
}
\details{A study was conducted in order to investigate a novel treatment
  in ankylosing spondylitis (Baeten et al., 2013). The primary endpoint
  related to \emph{treatment response}. 
  In order to formulate an informative prior distribution for the
  response rate to be expected in the control group of the new study, a
  systematic review of previous studies was consulted (McLeod et al.,
  2007), and, after a meta-analysis of the estimated response
  probabilities, the predictive distribution for the new study's
  response probability was derived. The predictive distribution here
  constitutes the \emph{meta-analytic-predictive (MAP) prior}
  distribution (Schmidli et al., 2014). The data set contains the
  relevant data from the eight \dQuote{historical} studies' placebo
  groups.

  Note that the original analysis (Baeten et al., 2013) involved a
  binomial model, and the resulting posterior predictive distribution
  was eventually approximated by a mixture of beta distributions.
}
\source{
  D. Baeten et al.
  Anti-interleukin-17A monoclonal antibody secukinumab in treatment of
  ankylosing spondylitis: a randomised, double-blind, placebo-controlled
  trial.
  \emph{The Lancet}, \bold{382}(9906):1705-1713, 2013.
  \doi{10.1016/S0140-6736(13)61134-4}.
}
\seealso{
  \code{\link{uisd}}, \code{\link{ess}}.
}
\references{
  C. McLeod et al.
  Adalimumab, etanercept, and infliximab for the treatment of ankylosing
  spondylitis: a systematic review and economic evaluation.
  \emph{Health Technology Assessment}, \bold{11}(28), 2007.
  \doi{10.3310/hta11280}.

  H. Schmidli, S. Gsteiger, S. Roychoudhury, A. O'Hagan,
  D. Spiegelhalter, B. Neuenschwander.
  Robust meta-analytic-predictive priors in clinical trials with
  historical control information.
  \emph{Biometrics}, \bold{70}(4):1023-1032, 2014.
  \doi{10.1111/biom.12242}.
}
\examples{
# load data:
data("BaetenEtAl2013")

# show data:
BaetenEtAl2013

\dontrun{
# compute effect sizes (logarithmic odds) from the count data:
as <- escalc(xi=events, ni=total, slab=study,
             measure="PLO", data=BaetenEtAl2013)

# compute the unit information standard deviation (UISD):
uisd(as)

# perform meta-analysis
# (using uniform priors for effect and heterogeneity):
bm <- bayesmeta(as)

# show results (log-odds):
forestplot(bm, xlab="log-odds", zero=NA)
# show results (odds):
forestplot(bm, exponentiate=TRUE, xlab="odds", zero=NA)

# show posterior predictive distribution --
# in terms of log-odds:
bm$summary[,"theta"]
# in terms of odds:
logodds <- exp(bm$summary[c(2,5,6), "theta"])
logodds
# in terms of probabilities:
(exp(logodds) / (exp(logodds) + 1))

# illustrate MAP prior density:
x <- seq(-3, 1, by=0.01)
plot(x, bm$dposterior(theta=x, predict=TRUE), type="l",
     xlab="log-odds (response)", ylab="posterior predictive density")
abline(h=0, v=0, col="grey")
}
}
\keyword{datasets}
