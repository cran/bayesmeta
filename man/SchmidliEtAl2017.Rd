\name{SchmidliEtAl2017}
\docType{data}
\alias{SchmidliEtAl2017}
\title{Historical variance example data}
\description{Estimates of endpoint variances from six studies.}
\usage{data("SchmidliEtAl2017")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{study} \tab \code{character} \tab study label \cr
    \bold{N}     \tab \code{integer}   \tab total sample size \cr
    \bold{stdev} \tab \code{numeric}   \tab standard deviation estimate \cr
    \bold{df}    \tab \code{integer}   \tab associated degrees of freedom
  }
}
\details{Schmidli \emph{et al.} (2017) investigated the use of
  information on an endpoint's variance from previous (\dQuote{historical})
  studies for the design and analysis of a new clinical trial. As an
  example application, the problem of designing a trial in \emph{wet
    age-related macular degeneration (AMD)} was considered. Trial
  design, and in particular considerations regarding the required sample
  size, hinge on the expected amount of variability in the primary
  endpoint (here: \emph{visual acuity}, which is measured on a scale
  ranging from 0 to 100 via an eye test chart).
  
  Historical data from six previous trials are available (Szabo \emph{et
  al.}; 2015), each trial providing an estimate \eqn{\hat{s}_i}{s[i]} of
  the endpoint's standard deviation along with the associated number of
  degrees of freedom \eqn{\nu_i}{nu[i]}. The standard deviations
  may then be modelled on the logarithmic scale, where the estimates and
  their associated standard errors are given by 
  \deqn{y_i=\log(\hat{s}_i) \quad \mbox{and} \quad 
  \sigma_i=\sqrt{\frac{1}{2\,\nu_i}}}{y[i] = log(s[i])  and  sigma[i] = sqrt(1/(2*nu[i]))}
  
  The \emph{unit information standard deviation (UISD)} for a logarithmic
  standard deviation then is at approximately
  \eqn{\frac{1}{\sqrt{2}}}{2^-0.5}. 
}
\source{
  H. Schmidli, B. Neuenschwander, T. Friede.
  Meta-analytic-predictive use of historical variance data
  for the design and analysis of clinical trials.
  \emph{Computational Statistics and Data Analysis}, \bold{113}:100-110, 2017.
  \doi{10.1016/j.csda.2016.08.007}.
}
\seealso{
  \code{\link{uisd}}, \code{\link{ess}}.
}
\references{
  S.M. Szabo, M. Hedegaard, K. Chan, K. Thorlund, R. Christensen,
  H. Vorum, J.P. Jansen.
  Ranibizumab vs. aflibercept for wet age-related macular degeneration:
  network meta-analysis to understand the value of reduced frequency dosing.
  \emph{Current Medical Research and Opinion}, \bold{31}(11):2031-2042, 2015.
  \doi{10.1185/03007995.2015.1084909}.
}
\examples{
# load data:
data("SchmidliEtAl2017")

# show data:
SchmidliEtAl2017

\dontrun{
# derive log-SDs and their standard errors:
dat <- cbind(SchmidliEtAl2017,
             logstdev    = log(SchmidliEtAl2017$stdev),
             logstdev.se = sqrt(0.5/SchmidliEtAl2017$df))
dat

# alternatively, use "metafor::escalc()" function:
es <- escalc(measure="SDLN",
             yi=log(stdev), vi=0.5/df, ni=N,
             slab=study, data=SchmidliEtAl2017)
es

# perform meta-analysis of log-stdevs:
bm <- bayesmeta(y=dat$logstdev,
                sigma=dat$logstdev.se,
                label=dat$study,
                tau.prior=function(t){dhalfnormal(t, scale=sqrt(2)/4)})

# or, alternatively:
bm <- bayesmeta(es,
                tau.prior=function(t){dhalfnormal(t, scale=sqrt(2)/4)})

# draw forest plot (see Fig.1):
forestplot(bm, zero=NA,
           xlab="log standard deviation")

# show heterogeneity posterior:
plot(bm, which=4, prior=TRUE)

# posterior of log-stdevs, heterogeneity,
# and predictive distribution:
bm$summary

# prediction (standard deviations):
exp(bm$summary[c(2,5,6),"theta"])
exp(bm$qposterior(theta=c(0.025, 0.25, 0.50, 0.75, 0.975), predict=TRUE))

# compute required sample size (per arm):
power.t.test(n=NULL, delta=8, sd=10.9, power=0.8)
power.t.test(n=NULL, delta=8, sd=14.0, power=0.8)

# check UISD:
uisd(es, indiv=TRUE)
uisd(es)
1 / sqrt(2)

# compute predictive distribution's ESS:
ess(bm, uisd=1/sqrt(2))
# actual total sample size:
sum(dat$N)

# illustrate predictive distribution
# on standard-deviation-scale (Fig.2):
x <- seq(from=5, to=20, length=200)
plot(x, (1/x) * bm$dposterior(theta=log(x), predict=TRUE), type="l",
     xlab=expression("predicted standard deviation "*sigma[k+1]),
     ylab="predictive density")
abline(h=0, col="grey")
}
}
\keyword{datasets}
