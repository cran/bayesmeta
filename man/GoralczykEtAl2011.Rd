\name{GoralczykEtAl2011}
\docType{data}
\alias{GoralczykEtAl2011}
\title{Liver transplant example data}
\description{Numbers of cases (transplant patients) and events (acute rejections,
  steroid resistant rejections, and deaths) in experimental and control groups of
  19 studies.}
\usage{data("GoralczykEtAl2011")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{publication}        \tab \code{character} \tab publication identifier (first author and publication year) \cr
    \bold{year}               \tab \code{numeric}   \tab publication year \cr
    \bold{randomized}         \tab \code{factor}    \tab randomization status (yes / no / not stated) \cr
    \bold{control.type}       \tab \code{factor}    \tab type of control group (\sQuote{concurrent} or \sQuote{historical})\cr
    \bold{comparison}         \tab \code{factor}    \tab type of comparison (\sQuote{IL-2RA only}, \sQuote{delayed CNI}, or \sQuote{no/low steroids}) \cr
    \bold{IL2RA}              \tab \code{factor}    \tab type of interleukin-2 receptor antagonist (IL-2RA) (\sQuote{basiliximab} or \sQuote{daclizumab}) \cr
    \bold{CNI}                \tab \code{factor}    \tab type of calcineurin inhibitor (CNI) (\sQuote{tacrolimus} or \sQuote{cyclosporine A}) \cr
    \bold{MMF}                \tab \code{factor}    \tab use of mycofenolate mofetil (MMF) (y/n) \cr
    \bold{followup}           \tab \code{numeric}   \tab follow-up time in months \cr
    \bold{treat.AR.events}    \tab \code{numeric}   \tab number of AR events in experimental group \cr
    \bold{treat.SRR.events}   \tab \code{numeric}   \tab number of SRR events in experimental group \cr
    \bold{treat.deaths}       \tab \code{numeric}   \tab number of deaths in experimental group \cr
    \bold{treat.total}        \tab \code{numeric}   \tab number of cases in experimental group \cr
    \bold{control.AR.events}  \tab \code{numeric}   \tab number of AR events in control group \cr
    \bold{control.SRR.events} \tab \code{numeric}   \tab number of SRR events in control group \cr
    \bold{control.deaths}     \tab \code{numeric}   \tab number of deaths in control group \cr
    \bold{control.total}      \tab \code{numeric}   \tab number of cases in control group \cr
  }
}
\details{A systematic literature review investigated the evidence on the
  effect of Interleukin-2 receptor antagonists (IL-2RA) and resulted in
  19 controlled studies reporting acute rejection (AR) and
  steroid-resistant rejection (SRR) rates as well as mortality in adult
  liver transplant recipients.
}
\source{A.D. Goralczyk, N. Hauke, N. Bari, T.Y. Tsui, T. Lorf, A. Obed.
  \href{https://doi.org/10.1002/hep.24385}{Interleukin-2 receptor
  antagonists for liver transplant recipients:
  A systematic review and meta-analysis of controlled studies}.
  \emph{Hepatology}, \bold{54}(2):541-554, 2011.
}
\seealso{
  \code{\link{CrinsEtAl2014}}.
}
\examples{
data("GoralczykEtAl2011")
\dontrun{
# compute effect sizes (log odds ratios) from count data
# (using "metafor" package's "escalc()" function):
require("metafor")
goralczyk.es <- escalc(measure="OR",
                       ai=exp.AR.events,  n1i=exp.total,
                       ci=cont.AR.events, n2i=cont.total,
                       slab=publication, data=GoralczykEtAl2011)
print(goralczyk.es[,c(1,10,12,13,15,16,17)])

# analyze using weakly informative half-Cauchy prior for heterogeneity:
goralczyk.ma <- bayesmeta(goralczyk.es, tau.prior=function(t){dhalfcauchy(t,scale=1)})

# show summary:
print(goralczyk.ma)

# show forest plot:
forestplot(goralczyk.ma)
}
}
\keyword{datasets}
