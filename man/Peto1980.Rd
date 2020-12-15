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
  \href{https://doi.org/10.1002/sim.650}{A comparison of statistical
    methods for meta-analysis}.
  \emph{Statistics in Medicine}, \bold{20}(6):825-840, 2001.
}
\references{R. Peto.
  \href{https://doi.org/10.1016/S0140-6736(80)91626-8}{Aspirin after
    myocardial infarction}. 
  \emph{The Lancet}, \bold{315}(8179):1172-1173, 1980.
  
  P.C. Elwood, A.L. Cochrane, M.L.Burr, P.M. Sweetnam, G. Williams,
  E. Welsby, S.J. Hughes, R. Renton.
  \href{https://doi.org/10.1136/bmj.1.5905.436}{A randomized controlled
  trial of acetyl salicylic acid in the secondary prevention of
  mortality from myocardial infarction}. 
  \emph{British Medical Journal}, \bold{1}(5905):436-440, 1974.

  The Coronary Drug Project Research Group.
  \href{https://doi.org/10.1016/0021-9681(76)90020-5}{Aspirin in coronary
  heart disease}.
  \emph{Journal of Chronic Diseases}, \bold{29}(10):625-642, 1976.

  K. Breddin, D. Loew, K. Lechner, K. Ueberla, E. Walter.
  \href{https://doi.org/10.1159/000214375}{Secondary prevention of
  myocardial infarction: a comparison of acetylsalicylic acid, placebo
  and phenprocoumon}.
  \emph{Haemostasis}, \bold{9}(6):325-344, 1980.

  P.C. Elwood, P.M. Sweetnam.
  \href{https://doi.org/10.1016/S0140-6736(79)92808-3}{Aspirin and
  secondary mortality after myocardial infarction}.
  \emph{The Lancet}, \bold{314}(8156):1313-1315, 1979.

  Aspirin Myocardial Infarction Study Research Group.
  \href{https://doi.org/10.1001/jama.1980.03300330019023}{A randomized,
  controlled trial of aspirin in persons recovered from myocardial
  infarction}.
  \emph{Journal of the American Medical Association},
  \bold{243}(7):661-669, 1980.

  The Persantine-Aspirin Reinfarction Study Research Group.
  \href{https://doi.org/10.1161/01.CIR.62.3.449}{Persantine and aspirin
    in coronary heart disease}.
  \emph{Circulation}, \bold{62}(3):449-461, 1980.  
}
\examples{
data("Peto1980")
\dontrun{
# compute effect sizes (log odds ratios) from count data
# (using "metafor" package's "escalc()" function):
require("metafor")
peto.es <- escalc(measure="OR",
                  ai=treat.events,   n1i=treat.cases,
                  ci=control.events, n2i=control.cases,
                  slab=publication, data=Peto1980)
print(peto.es)

# check sensitivity to different prior choices:
peto.ma01 <- bayesmeta(peto.es)
peto.ma02 <- bayesmeta(peto.es, tau.prior=function(t){dhalfnormal(t, scale=1)})
 
par(mfrow=c(2,1))
  plot(peto.ma01, which=4, prior=TRUE, taulim=c(0,1), main="uniform prior")
  plot(peto.ma02, which=4, prior=TRUE, taulim=c(0,1), main="half-normal prior")
par(mfrow=c(1,1))

# compare heterogeneity (tau) estimates:
print(rbind("uniform"    =peto.ma01$summary[,"tau"],
            "half-normal"=peto.ma02$summary[,"tau"]))

# compare effect (mu) estimates:
print(rbind("uniform"    =peto.ma01$summary[,"mu"],
            "half-normal"=peto.ma02$summary[,"mu"]))

summary(peto.ma02)
forestplot(peto.ma02)
plot(peto.ma02)
}
}
\keyword{datasets}