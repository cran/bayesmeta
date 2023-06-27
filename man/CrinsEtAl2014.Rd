\name{CrinsEtAl2014}
\docType{data}
\alias{CrinsEtAl2014}
\title{Pediatric liver transplant example data}
\description{Numbers of cases (transplant patients) and events (acute rejections,
  steroid resistant rejections, PTLDs, and deaths) in experimental and control groups of
  six studies.}
\usage{data("CrinsEtAl2014")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{publication}         \tab \code{character} \tab publication identifier (first author and publication year) \cr
    \bold{year}                \tab \code{numeric}   \tab publication year \cr
    \bold{randomized}          \tab \code{factor}    \tab randomization status (y/n) \cr
    \bold{control.type}        \tab \code{factor}    \tab type of control group (\sQuote{concurrent} or \sQuote{historical})\cr
    \bold{comparison}          \tab \code{factor}    \tab type of comparison (\sQuote{IL-2RA only}, \sQuote{delayed CNI}, or \sQuote{no/low steroids}) \cr
    \bold{IL2RA}               \tab \code{factor}    \tab type of interleukin-2 receptor antagonist (IL-2RA) (\sQuote{basiliximab} or \sQuote{daclizumab}) \cr
    \bold{CNI}                 \tab \code{factor}    \tab type of calcineurin inhibitor (CNI) (\sQuote{tacrolimus} or \sQuote{cyclosporine A}) \cr
    \bold{MMF}                 \tab \code{factor}    \tab use of mycofenolate mofetil (MMF) (y/n) \cr
    \bold{followup}            \tab \code{numeric}   \tab follow-up time in months \cr
    \bold{treat.AR.events}     \tab \code{numeric}   \tab number of AR events in experimental group \cr
    \bold{treat.SRR.events}    \tab \code{numeric}   \tab number of SRR events in experimental group \cr
    \bold{treat.PTLD.events}   \tab \code{numeric}   \tab number of PTLD events in experimental group \cr
    \bold{treat.deaths}        \tab \code{numeric}   \tab number of deaths in experimental group \cr
    \bold{treat.total}         \tab \code{numeric}   \tab number of cases in experimental group \cr
    \bold{control.AR.events}   \tab \code{numeric}   \tab number of AR events in control group \cr
    \bold{control.SRR.events}  \tab \code{numeric}   \tab number of SRR events in control group \cr
    \bold{control.PTLD.events} \tab \code{numeric}   \tab number of PTLD events in control group \cr
    \bold{control.deaths}      \tab \code{numeric}   \tab number of deaths in control group \cr
    \bold{control.total}       \tab \code{numeric}   \tab number of cases in control group \cr
  }
}
\details{A systematic literature review investigated the evidence on the
  effect of Interleukin-2 receptor antagonists (IL-2RA) and resulted in
  six controlled studies reporting acute rejection (AR),
  steroid-resistant rejection (SRR) and post-transplant
  lymphoproliferative disorder (PTLD) rates as well as mortality in
  pediatric liver transplant recipients.
}
\source{N.D. Crins, C. Roever, A.D. Goralczyk, T. Friede.
  Interleukin-2 receptor antagonists for pediatric liver transplant
  recipients: A systematic review and meta-analysis of controlled
  studies.
  \emph{Pediatric Transplantation}, \bold{18}(8):839-850, 2014.
  \doi{10.1111/petr.12362}.
}
\seealso{
  \code{\link{GoralczykEtAl2011}}.
}
\references{
  C. Roever.
  Bayesian random-effects meta-analysis using the bayesmeta R package.
  \emph{Journal of Statistical Software}, \bold{93}(6):1-51, 2020.
  \doi{10.18637/jss.v093.i06}.
  
  C. Roever, T. Friede.
  Using the bayesmeta R package for Bayesian random-effects meta-regression.
  \emph{Computer Methods and Programs in Biomedicine},
  \bold{299}:107303, 2023.
  \doi{10.1016/j.cmpb.2022.107303}.

  T.G. Heffron et al.
  Pediatric liver transplantation with daclizumab induction therapy.
  \emph{Transplantation}, \bold{75}(12):2040-2043, 2003.
  \doi{10.1097/01.TP.0000065740.69296.DA}.
  
  N.E.M. Gibelli et al.
  Basiliximab-chimeric anti-IL2-R monoclonal antibody in pediatric liver
  transplantation: comparative study.
  \emph{Transplantation Proceedings}, \bold{36}(4):956-957, 2004.
  \doi{10.1016/j.transproceed.2004.04.070}.

  S. Schuller et al.
  Daclizumab induction therapy associated with tacrolimus-MMF has better
  outcome compared with tacrolimus-MMF alone in pediatric living donor
  liver transplantation.
  \emph{Transplantation Proceedings}, \bold{37}(2):1151-1152, 2005.
  \doi{10.1016/j.transproceed.2005.01.023}.

  R. Ganschow et al.
  Long-term results of basiliximab induction immunosuppression in
  pediatric liver transplant recipients.
  \emph{Pediatric Transplantation}, \bold{9}(6):741-745, 2005.
  \doi{10.1111/j.1399-3046.2005.00371.x}.

  M. Spada et al.
  Randomized trial of basiliximab induction versus steroid therapy in
  pediatric liver allograft recipients under tacrolimus immunosuppression.
  \emph{American Journal of Transplantation}, \bold{6}(8):1913-1921, 2006.
  \doi{10.1111/j.1600-6143.2006.01406.x}.
  
  J.M. Gras et al.
  Steroid-free, tacrolimus-basiliximab immunosuppression in pediatric
  liver transplantation: Clinical and pharmacoeconomic study in 50 children.
  \emph{Liver Transplantation}, \bold{14}(4):469-477, 2008.
  \doi{10.1002/lt.21397}.
}
\examples{
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
