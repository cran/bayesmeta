\name{Peto1980}
\docType{data}
\alias{Peto1980}
\title{Aspirin after myocardial infarction example data}
\description{Numbers of cases (patients) and events (deaths) in
  treatment and control groups of six studies.}
\usage{data("Peto1980")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{publication}    \tab \code{character} \tab publication reference \cr
    \bold{study}          \tab \code{character} \tab study acronym or abbreviation \cr
    \bold{start}, \bold{end} \tab \code{integer}   \tab duration of study (calendar years) \cr
    \bold{age}            \tab \code{numeric}   \tab mean patient age (years) \cr
    \bold{dose}           \tab \code{numeric}   \tab total daily dose (mg) \cr
    \bold{followup}       \tab \code{numeric}   \tab follow-up duration (months) \cr
    \bold{treat.cases}    \tab \code{integer}   \tab number of cases in treatment group \cr
    \bold{treat.events}   \tab \code{integer}   \tab number of events in treatment group \cr
    \bold{control.cases}  \tab \code{integer}   \tab number of cases in control group \cr
    \bold{control.events} \tab \code{integer}   \tab number of events in control group \cr
  }
}
\details{
  Peto (1980) investigated mortality data from six randomized,
  placebo-controlled clinical trials of aspirin, involving a total of
  10,703 post-myocardial infarction patients. Canner (1987) later investigated
  potential heterogeneity between study characteristics as well as their
  reported estimates. The included studies' abbreviations are:
  \tabular{ll}{
    UK-1   \tab  first United Kingdom trial \cr
    CDPA   \tab  Coronary Drug Project Aspirin trial \cr
    GAMS   \tab  German-Austrian Multicentre Study\cr
    UK-2   \tab  second United Kingdom trial \cr
    PARIS  \tab  Persantine-Aspirin Reinfarction Study \cr
    AMIS   \tab  Aspirin Myocardial Infarction Study
  }
}
\source{P.L. Canner. An overview of six clinical trials of aspirin
  in coronary heart disease.
  \emph{Statistics in Medicine}, \bold{6}(3):255-263, 1987.
  \doi{10.1002/sim.4780060310}
}
\references{R. Peto.
  Aspirin after myocardial infarction. 
  \emph{The Lancet}, \bold{315}(8179):1172-1173, 1980.
  \doi{10.1016/S0140-6736(80)91626-8}.
  
  P.C. Elwood, A.L. Cochrane, M.L.Burr, P.M. Sweetnam, G. Williams,
  E. Welsby, S.J. Hughes, R. Renton.
  A randomized controlled trial of acetyl salicylic acid in the
  secondary prevention of mortality from myocardial infarction. 
  \emph{British Medical Journal}, \bold{1}(5905):436-440, 1974.
  \doi{10.1136/bmj.1.5905.436}.

  The Coronary Drug Project Research Group.
  Aspirin in coronary heart disease.
  \emph{Journal of Chronic Diseases}, \bold{29}(10):625-642, 1976.
  \doi{10.1016/0021-9681(76)90020-5}.

  K. Breddin, D. Loew, K. Lechner, K. Ueberla, E. Walter.
  Secondary prevention of myocardial infarction: a comparison of
  acetylsalicylic acid, placebo and phenprocoumon.
  \emph{Haemostasis}, \bold{9}(6):325-344, 1980.
  \doi{10.1159/000214375}.

  P.C. Elwood, P.M. Sweetnam.
  Aspirin and secondary mortality after myocardial infarction.
  \emph{The Lancet}, \bold{314}(8156):1313-1315, 1979.
  \doi{10.1016/S0140-6736(79)92808-3}.

  Aspirin Myocardial Infarction Study Research Group.
  A randomized, controlled trial of aspirin in persons recovered from
  myocardial infarction.
  \emph{Journal of the American Medical Association},
  \bold{243}(7):661-669, 1980.
  \doi{10.1001/jama.1980.03300330019023}.

  The Persantine-Aspirin Reinfarction Study Research Group.
  Persantine and aspirin in coronary heart disease.
  \emph{Circulation}, \bold{62}(3):449-461, 1980.
  \doi{10.1161/01.CIR.62.3.449}.
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