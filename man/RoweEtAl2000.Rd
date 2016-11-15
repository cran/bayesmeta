\name{RoweEtAl2000}
\docType{data}
\alias{RoweEtAl2000}
\title{Magnesium sulfate example data}
\description{Numbers of hospitalisation events in asthma patients treated with intravenous magnesium sulfate.}
\usage{data("RoweEtAl2000")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{study}          \tab \code{character} \tab publication identifier (first author and publication year) \cr
    \bold{age.range}      \tab \code{character} \tab age range of included patients \cr
    \bold{events.MgSO4}   \tab \code{numeric}   \tab number of patients admitted to hospital in experimental group \cr
    \bold{total.MgSO4}    \tab \code{numeric}   \tab total number of patients in experimental group \cr
    \bold{events.placebo} \tab \code{numeric}   \tab number of patients admitted to hospital in control group \cr
    \bold{total.placebo}  \tab \code{numeric}   \tab total number of patients in control group \cr
  }
}
\details{
  A systematic literature search resulted in four randomized,
  placebo-controlled studies in patients with severe acute asthma 
  and treated with magnesium sulfate (\eqn{\mathrm{MgSO}_4}{MgSO4}) versus placebo
  that reported data on hopital admissions.
  Two of the studies (Skobeloff (1989) and Bloch (1985)) were based on
  adults, while the other two (Ciarallo (1997) and Devi (1997))
  investigated the effect in children.
}
\source{
  B.H. Rowe, J. Bretzlaff, C. Bourdon, S. Blitz, C.A. Camargo.
  \href{http://dx.doi.org/10.1002/14651858.CD001490}{Magnesium 
  sulfate for treating exacerbations of acute asthma 
  in the emergency department}.
  \emph{Cochrane Database of Systematic Reviews}, 
  \bold{1}:CD001490, 2000.
}
\references{
  E.M. Skobeloff, W. H. Spivey, R.M. McNamara, L. Greenspon.
  \href{http://dx.doi.org/10.1001/jama.1989.03430090072036}{Intravenous
  magnesium sulfate for the treatment of acute asthma in the emergency
  department}. \emph{JAMA}, \bold{262}(9):1210-1213, 1989.

  H. Bloch, R. Silverman, N. Mancherje, S. Grant, L. Jagminas,
  S.M. Scharf.
  \href{http://dx.doi.org/10.1378/chest.107.6.1576}{Intravenous
  magnesium sulfate as an adjunct in the treatment of acute asthma}.
  \emph{Chest}, \bold{107}(6):1576-1581, 1995.

  L. Ciarallo, A.H. Sauer, M.W. Shannon.
  \href{http://dx.doi.org/10.1016/S0022-3476(96)70023-9}{Intravenous
  magnesium therapy for moderate to severe pediatric asthma: Results
  of a randomized, placebo-controlled trial}
  \emph{The Journal of Pediatrics}, \bold{129}(6):809-814, 1996.

  P.R. Devi, L. Kumar, S.C. Singhi, R. Prasad, M. Singh.
  Intravenous magnesium sulfate in acute severe asthma 
  not responding to conventional therapy.
  \emph{Indian Pediatrics}, \bold{34}(5):389-397, 1997.
}
\examples{
# load data:
data("RoweEtAl2000")

# show data:
print(RoweEtAl2000)

\dontrun{
# compute effect sizes
# (log-ORs, adult studies only):
if (require("metafor")) {
  es.adults   <- escalc(measure="OR",
                        ai=events.MgSO4,  n1i=total.MgSO4,
                        ci=events.placebo, n2i=total.placebo,
                        slab=study, data=RoweEtAl2000[1:2,])
} else print("Sorry, 'metafor' package not installed!")

# show effect sizes:
print(es.adults)

# analyze data:
ma.adults <- bayesmeta(es.adults,
                       mu.prior.mean=0, mu.prior.sd=3,
                       tau.prior=function(t){dhalfcauchy(t,scale=0.5)})

# show posterior distribution of effect (log-OR):
plot(ma.adults, which=3, prior=TRUE)
abline(v=0, lty="dotted")

forestplot(ma.adults)
}
}
\keyword{datasets}
