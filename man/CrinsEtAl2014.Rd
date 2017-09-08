\name{CrinsEtAl2014}
\docType{data}
\alias{CrinsEtAl2014}
\title{Pediatric liver transplant example data}
\description{Numbers of cases (transplants) and events (acute rejections
  or steroid resistant rejections) in experimental and control groups of
  six studies.}
\usage{data("CrinsEtAl2014")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{publication}        \tab \code{character} \tab publication identifier (first author and publication year) \cr
    \bold{year}               \tab \code{numeric}   \tab publication year \cr
    \bold{randomized}         \tab \code{factor}    \tab randomization status (y/n) \cr
    \bold{control.type}       \tab \code{factor}    \tab type of control group (\sQuote{concurrent} or \sQuote{historical})\cr
    \bold{comparison}         \tab \code{factor}    \tab type of comparison (\sQuote{IL-2RA only}, \sQuote{delayed CNI}, or \sQuote{no/low steroids}) \cr
    \bold{IL2RA}              \tab \code{factor}    \tab type of interleukin-2 receptor antagonist (IL-2RA) (\sQuote{basiliximab} or \sQuote{daclizumab}) \cr
    \bold{CNI}                \tab \code{factor}    \tab type of calcineurin inhibitor (CNI) (\sQuote{tracrolimus} or \sQuote{cyclosporine A}) \cr
    \bold{MMF}                \tab \code{factor}    \tab use of mycofenolate mofetil (MMF) (y/n) \cr
    \bold{followup}           \tab \code{numeric}   \tab follow-up time in months \cr
    \bold{treat.AR.events}    \tab \code{numeric}   \tab number of AR events in experimental group \cr
    \bold{treat.SRR.events}   \tab \code{numeric}   \tab number of SRR events in experimental group \cr
    \bold{treat.total}        \tab \code{numeric}   \tab number of cases in experimental group \cr
    \bold{control.AR.events}  \tab \code{numeric}   \tab number of AR events in control group \cr
    \bold{control.SRR.events} \tab \code{numeric}   \tab number of SRR events in control group \cr
    \bold{control.total}      \tab \code{numeric}   \tab number of cases in control group \cr
  }
}
\details{A systematic literature review investigated the evidence on
  the effect of Interleukin-2 receptor antagonists (IL-2RA) and resulted
  in six controlled studies reporting acute rejection (AR) and
  steroid-resistant rejection (SRR) rates in pediatric liver transplant
  recipients.
}
\source{N.D. Crins, C. Roever, A.D. Goralczyk, T. Friede.
  \href{http://doi.org/10.1111/petr.12362}{Interleukin-2 receptor
  antagonists for pediatric liver transplant recipients:
  A systematic review and meta-analysis of controlled studies}.
  \emph{Pediatric Transplantation}, \bold{18}(8):839-850, 2014.
}
\references{T.G. Heffron et al.
  Pediatric liver transplantation with daclizumab induction therapy.
  \emph{Transplantation}, \bold{75}(12):2040-2043, 2003.
  
  N.E.M. Gibelli et al.
  \href{http://doi.org/10.1016/j.transproceed.2004.04.070}{Basiliximab-chimeric
    anti-IL2-R monoclonal antibody in pediatric liver transplantation:
    comparative study}.
  \emph{Transplantation Proceedings}, \bold{36}(4):956-957, 2004.

  S. Schuller et al.
  \href{http://doi.org/10.1016/j.transproceed.2005.01.023}{Daclizumab 
    induction therapy associated with tacrolimus-MMF has better outcome
  compared with tacrolimus-MMF alone in pediatric living donor liver
  transplantation}.
  \emph{Transplantation Proceedings}, \bold{37}(2):1151-1152, 2005.

  R. Ganschow et al.
  \href{http://doi.org/10.1111/j.1399-3046.2005.00371.x}{Long-term
  results of basiliximab induction immunosuppression in pediatric liver
  transplant recipients}.
  \emph{Pediatric Transplantation}, \bold{9}(6):741-745, 2005.

  M. Spada et al.
  \href{http://doi.org/10.1111/j.1600-6143.2006.01406.x}{Randomized
  trial of basiliximab induction versus steroid therapy in pediatric
  liver allograft recipients under tacrolimus
  immunosuppression}.
  \emph{American Journal of Transplantation}, \bold{6}(8):1913-1921, 2006.

  J.M. Gras et al.
  \href{http://doi.org/10.1002/lt.21397}{Steroid-free,
  tacrolimus-basiliximab immunosuppression in pediatric liver
  transplantation: Clinical and pharmacoeconomic study in 50 children}.
  \emph{Liver Transplantation}, \bold{14}(4):469-477, 2008.}
\examples{
data("CrinsEtAl2014")
\dontrun{
# compute effect sizes (log odds ratios) from count data:
if (require("metafor")) {
  crins.es <- escalc(measure="OR",
                     ai=exp.AR.events,  n1i=exp.total,
                     ci=cont.AR.events, n2i=cont.total,
                     slab=publication, data=CrinsEtAl2014)
} else print("Sorry, 'metafor' package not installed!")
print(crins.es)

# analyze using weakly informative half-Cauchy prior for heterogeneity:
crins.ma <- bayesmeta(crins.es, tau.prior=function(t){dhalfcauchy(t,scale=1)})

# show results:
print(crins.ma)
forestplot(crins.ma)
plot(crins.ma)

# show heterogeneity posterior along with prior:
plot(crins.ma, which=4, prior=TRUE)

# perform meta analysis using 2 randomized studies only
# but use 4 non-randomized studies to inform heterogeneity prior:
crins.nrand <- bayesmeta(crins.es[crins.es$randomized=="no",],
                         tau.prior=function(t){dhalfcauchy(t,scale=1)})
crins.rand  <- bayesmeta(crins.es[crins.es$randomized=="yes",],
                         tau.prior=function(t){crins.nrand$dposterior(tau=t)})
plot(crins.nrand, which=4, prior=TRUE,
     main="non-randomized posterior = randomized prior")
plot(crins.rand, which=4, prior=TRUE, main="randomized posterior")
plot(crins.rand, which=1)
}
}
\keyword{datasets}
