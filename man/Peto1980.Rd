\name{Peto1980}
\docType{data}
\alias{Peto1980}
\title{Aspirin after myocardial infarction example data}
\description{Numbers of cases (patients) and events (deaths) in
  treatment and control groups of six studies.}
\usage{data("Peto1980")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{publication}    \tab \code{character} \tab publication identifier \cr
    \bold{treat.cases}    \tab \code{numeric}   \tab number of cases in treatment group \cr
    \bold{treat.events}   \tab \code{numeric}   \tab number of events in treatment group \cr
    \bold{control.cases}  \tab \code{numeric}   \tab number of cases in control group \cr
    \bold{control.events} \tab \code{numeric}   \tab number of events in control group \cr
  }
}
\details{
  Quoting from Brockwell and Gordon (2001):
  \dQuote{The collection consists of six studies, each examining the
    effect of aspirin after myocardial infarction. In each study the
    number of patients who died after having been given either aspirin
    or a control drug is recorded.}
}
\source{S.E. Brockwell, I.R. Gordon.
  \href{http://dx.doi.org/10.1002/sim.650}{A comparison of statistical
    methods for meta-analysis}.
  \emph{Statistics in Medicine}, \bold{20}(6):825-840, 2001.
}
\references{R. Peto.
  \href{http://dx.doi.org/10.1016/S0140-6736(80)91626-8}{Aspirin after
    myocardial infarction}. 
  \emph{The Lancet}, \bold{315}(8179):1172-1173, 1980.
}
\examples{
\dontrun{
data("Peto1980")
if (require("metafor")) {
  peto.es <- escalc(measure="OR",
                    ai=treat.events,   n1i=treat.cases,
                    ci=control.events, n2i=control.cases,
                    slab=publication, data=Peto1980)
} else print("Sorry, 'metafor' package not installed!")
print(peto.es)

# check sensitivity to different prior choices:
peto.ma01 <- bayesmeta(peto.es)
peto.ma02 <- bayesmeta(peto.es, tau.prior=function(t){dhalfnormal(t, scale=1)})
 
par(mfrow=c(2,1))
  plot(peto.ma01, which=4, prior=TRUE, main="uniform prior")
  plot(peto.ma02, which=4, prior=TRUE, main="half-normal prior")
par(mfrow=c(1,1))

# compare heterogeneity (tau) estimates:
print(rbind("uniform"    =peto.ma01$summary[,"tau"],
            "half-normal"=peto.ma02$summary[,"tau"]))

# compare effect (mu) estimates:
print(rbind("uniform"    =peto.ma01$summary[,"mu"],
            "half-normal"=peto.ma02$summary[,"mu"]))

summary(peto.ma02)
plot(peto.ma02)
}
}
\keyword{datasets}