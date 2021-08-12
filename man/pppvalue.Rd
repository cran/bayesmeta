\name{pppvalue}
\alias{pppvalue}
\title{Posterior predictive \eqn{p}-values}
\description{
  Compute posterior or prior predictive \eqn{p}-values from a
  \code{\link{bayesmeta}} object.
}
\usage{
  pppvalue(x, parameter = "mu", value = 0.0,
           alternative = c("two.sided", "less", "greater"),
           statistic = "median",
           rejection.region,
           n = 10,
           prior = FALSE,
           quietly = FALSE,
           parallel, seed, ...)
}
\arguments{
  \item{x}{a \code{bayesmeta} object.}
  \item{parameter}{the parameter to be tested. May be the effect
    (\code{"mu"}), the heterogeneity (\code{"tau"}) or one of the
    study-specific (\eqn{\theta_i}{theta[i]}) parameters denoted by
    their label or their index.}
  \item{value}{the (null-) hypothesized value.}
  \item{alternative}{the type of alternative hypothesis.}
  \item{statistic}{the figure to be used a the \sQuote{test
    statistic}, or \sQuote{discrepancy variable}. May be chosen as \code{"t"},
    \code{"Q"} or \code{"cdf"}, or among the row
    names of the \code{bayesmeta} object's \sQuote{\code{...$summary}}
    element. \emph{Or} it may be specified as a \code{function}. For details, see below.}
  \item{rejection.region}{the test statistic's rejection region. May be
    one of \code{"upper.tail"}, \code{"lower.tail"} or
    \code{"two.tailed"}. If unspecified, it is set automatically based
    on the \sQuote{\code{alternative}} and \sQuote{\code{statistic}}
    parameters.}
  \item{n}{the number of Monte Carlo replications to be generated. The
    default value is \code{n=10}, but in practice a substantially larger
    value should be appropriate.}
  \item{prior}{a logical flag to request \emph{prior predictive} (instead
    of \emph{posterior predictive}) \eqn{p}-values. Prior predictive
    values are only available for hypotheses concerning the effect
    (\eqn{\mu}) and heterogeneity (\eqn{\tau}) parameters.}
  \item{quietly}{a logical flag to show (or suppress) output during
    computation; this may also speed up computations slightly.}
  \item{parallel}{the number of parallel processes to utilize. By
    default, if multiple (k) cores are detected, then k-1 parallel
    processes are used.}
  \item{seed}{(optional) an \code{integer} random seed value to
    generate reproducible results.}
  \item{...}{further parameters passed to \sQuote{\code{statistic}},
  \emph{if} the \sQuote{\code{statistic}} argument was specified as a
  \code{function}.}
}
\details{
  Posterior predictive \eqn{p}-values are Bayesian analogues to
  \sQuote{classical} \eqn{p}-values (Meng, 1994; Gelman, Meng and Stern,
  1996; Gelman et al., 2014). The \code{pppvalue()} function allows to
  compute these values for one- and two-sided hypotheses concerning the
  effect (\eqn{\mu}) or heterogeneity (\eqn{\tau}) parameter, or one of
  the study-specific effect parameters (\eqn{\theta_i}{theta[i]}) in a
  random-effects meta-analysis.

  \emph{Prior} predictive \eqn{p}-values have a
  similar interpretation, but they have a stronger dependence on the
  prior specification and are only available when the prior is proper;
  for a more detailed discussion, see Gelman, Meng and Stern (1996,
  Sec. 4).

  The function may also be used to only generate samples (\eqn{\tau},
  \eqn{\mu}, \eqn{\theta_i}{theta[i]}, \eqn{y_i}{y[i]}) without having to also
  derive a statistic or a \eqn{p}-value. In order to achieve that, the
  \sQuote{\code{statistic}} argument may be specified as
  \sQuote{\code{NA}}, and generated samples may be recovered from the
  \sQuote{\code{...$replicates}} output element.
  
  \subsection{\eqn{p}-values from Monte Carlo sampling}{The computation
  is based on Monte Carlo sampling and repeated analysis of re-generated
  data sets drawn from the parameters' (conditional) posterior
  predictive (or prior) distribution; so the \eqn{p}-value derivation is
  somewhat computationally expensive. The \eqn{p}-value eventually is
  computed based on how far in the tail area the actual data are (in
  terms of the realized \sQuote{test statistic} or \sQuote{discrepancy})
  relative to the Monte-Carlo-sampled distribution. Accuracy of the
  computed \eqn{p}-value hence strongly depends on the number of samples
  (as specified through the \sQuote{\code{n}} argument) that are
  generated. Also, (near-) zero \eqn{p}-values need to be interpreted
  with caution, and in relation to the used Monte Carlo sample size
  (\code{n}).}

  \subsection{\sQuote{Test}-statistics or \sQuote{discrepancy
  variables}}{The \sQuote{\code{statistic}} argument determines the statistic
  to be computed from the data as a measure of deviation from the null
  hypothesis. If specified as \code{"Q"}, then Cochran's \eqn{Q} statistic is
  computed; this is useful for testing for homogeneity (\eqn{\tau=0}). If specified as
  one of the row names of the \sQuote{\code{x$summary}} element, then,
  depending on the type of null hypothesis specified through the
  \sQuote{\code{parameter}} argument, the corresponding parameter's posterior
  quantity is used for the statistic. If specified as \code{"t"}, then a
  \eqn{t}-type statistic is computed (the difference between the
  corresponding parameter's posterior mean and its hypothesized value,
  divided by the posterior standard deviation). If specified as
  \code{"cdf"}, the parameter's marginal posterior cumulative
  distribution function evaluated a the hypothesized value
  (\sQuote{\code{value}}) is used.

  The \sQuote{\code{statistic}} argument may also be specified as an
  arbitrary \code{function} of the data (\eqn{y}). The \code{function}'s
  first argument then needs to be the data (\eqn{y}), additional
  arguments may be passed as arguments (\sQuote{\code{...}}) to the
  \sQuote{\code{pppvalue()}} function. See also the examples below.}

  \subsection{One- and two-sided hypotheses}{Specification of one- or
  two-sided hypotheses not only has implications for the determination
  of the \eqn{p}-value from the samples, but also for the sampling
  process itself. Parameter values are drawn from a subspace according
  to the null hypothesis, which for a two-sided test is a line, and for
  a one-sided test is a half-plane. This also implies that one- and
  two-sided \eqn{p}-values cannot simply be converted into one
  another.

  For example, when specifying
  \code{pppvalue(..., param="mu", val=0, alt="two.sided")},
  then first paramater values (\eqn{\tau}, \eqn{\mu}) are drawn from the
  conditional posterior distribution \eqn{p(\tau, \mu | y, \sigma,
  \mu=0)}{p(tau, mu | y, sigma, mu=0)}, and subsequently new data sets
  are generated based on the parameters. If a one-sided hypothesis is
  specified, e.g. via
  \code{pppvalue(..., param="mu", val=0, alt="less")},
  then parameters are drawn from \eqn{p(\tau, \mu | y,
  \sigma, \mu>0)}{p(tau, mu | y, sigma, mu>0)}.

  For a hypothesis concerning the individual effect parameters
  \eqn{\theta_i}{theta[i]}, conditions are imposed on the corresponding
  \eqn{\theta_i}{theta[i]}. For example, for a specification of
  \code{pppvalue(..., param=2, val=0, alt="less")}, the
  hypothesis concerns the \eqn{i}=2nd study's effect paramater
  \eqn{\theta_2}{theta[2]}. First a sample is generated from
  \eqn{p(\theta_2|y, \sigma, \theta_2 > 0)}{p(theta[2] | y, sigma,
  theta[2]>0)}. Then samples of \eqn{\mu} and \eqn{\tau} are generated
  by conditioning on the generated \eqn{\theta_2}{theta[2]} value, and
  data \eqn{y} are generated by conditioning on all three.

  Unless explicitly specified through the
  \sQuote{\code{rejection.region}} argument, the test statistic's
  \dQuote{rejection region} (the direction in which extreme \code{statistic}
  values indicate a departure from the null hypothesis) is set based on the
  \sQuote{\code{alternative}} and \sQuote{\code{statistic}}
  parameters. The eventually used setting can be checked in the output's
  \sQuote{\code{...$rejection.region}} component.}

  \subsection{Computation}{When aiming to compute a \eqn{p}-value, it is
  probably a good idea to first start with a smaller \sQuote{\code{n}}
  argument to get a rough idea of the \eqn{p}-value's order of magnitude
  as well as the computational speed, before going over to a larger,
  more realistic \code{n} value. The implementation is able to utilize
  multiple processors or cores via the \pkg{parallel} package; details
  may be specified via the \sQuote{\code{parallel}} argument.}

}
\value{
  A \code{list} of class \sQuote{\code{htest}} containing the following
  components:
  \item{statistic}{the \sQuote{test statistic} (or \sQuote{discrepancy})
    value based on the actual data.}
  \item{parameter}{the number (\code{n}) of Monte Carlo replications used.}
  \item{p.value}{the derived \eqn{p}-value.}
  \item{null.value}{the (null-) hypothesized parameter value.}
  \item{alternative}{the type of alternative hypothesis.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{the name of the underlying \code{bayesmeta} object.}
  \item{call}{an object of class \code{call} giving the
    function call that generated the \code{htest} object.}
  \item{rejection.region}{the test statistic's rejection region.}
  \item{replicates}{a \code{list} containing the replicated parameters
    (\eqn{\tau}, \eqn{\mu}, \eqn{\theta_i}{theta[i]}),
    data (\eqn{y_i}{y[i]}) and statistic, along with an indicator for those
    samples constituting the distribution's \sQuote{tail area}.}
  \item{computation.time}{The computation time (in seconds) used.}
}
\references{
  X.-L. Meng.
  Posterior predictive p-values.
  \emph{The Annals of Statistics}, \bold{22}(3):1142-1160, 1994.
  \doi{10.1214/aos/1176325622}.

  A. Gelman, X.-L. Meng, H. Stern.
  Posterior predictive assessment of model fitness
  via realized discrepancies.
  \emph{Statistica Sinica}, \bold{6}(4):733-760, 1996.
  
  A. Gelman, J.B. Carlin, H.S. Stern, D.B. Dunson, A. Vehtari,
  D.B. Rubin.
  \emph{Bayesian data analysis}.
  Chapman & Hall / CRC, Boca Raton, 2014.
  
  C. Roever.
  Bayesian random-effects meta-analysis using the bayesmeta R package.
  \emph{Journal of Statistical Software}, \bold{93}(6):1-51, 2020.
  \doi{10.18637/jss.v093.i06}.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{bayesmeta}}, \code{\link[stats]{prop.test}}.
}
\examples{
\dontrun{
# perform a meta analysis;
# load data:
data("CrinsEtAl2014")
# compute effect sizes (log odds ratios) from count data
# (using "metafor" package's "escalc()" function):
require("metafor")
crins.srr <- escalc(measure="OR",
                    ai=exp.SRR.events,  n1i=exp.total,
                    ci=cont.SRR.events, n2i=cont.total,
                    slab=publication, data=CrinsEtAl2014, subset=c(1,4,6))
# analyze:
bma <- bayesmeta(crins.srr, mu.prior.mean=0, mu.prior.sd=4,
                 tau.prior=function(t){dhalfnormal(t, scale=0.5)})

# compute a 2-sided p-value for the effect (mu) parameter
# (note: this may take a while!):
p <- pppvalue(bma, parameter="mu", value=0, n=100)

# show result:
print(p)

# show the test statistic's distribution
# along with its actualized value:
plot(ecdf(p$replicates$statistic[,1]),
     xlim=range(c(p$statistic, p$replicates$statistic[,1])))
abline(v=p$statistic, col="red")

# show the parameter values
# drawn from the (conditional) posterior distribution:
plot(bma, which=2)
abline(h=p$null.value)                                # (the null-hypothesized mu value)
points(p$replicates$tau, p$replicates$mu, col="cyan") # (the samples)

######################################################################
#  Among the 3 studies, only the first (Heffron, 2003) was randomized.
#  One might wonder about this particular study's effect (theta[1])
#  in the light of the additional evidence and compute a one-sided
#  p-value:

p <- pppvalue(bma, parameter="Heffron", value=0, n=100, alternative="less")
print(p)

######################################################################
#  One may also define one's own 'test' statistic to be used.
#  For example, one could utilize the Bayes factor to generate
#  a p-value for the homogeneity (tau=0) hypothesis:

BF <- function(y, sigma)
{
  bm <- bayesmeta(y=y, sigma=sigma,
                  mu.prior.mean=0, mu.prior.sd=4,
                  tau.prior=function(t){dhalfnormal(t, scale=0.5)},
                  interval.type="central")
  # (central intervals are faster to compute;
  #  interval type otherwise is not relevant here)
  return(bm$bayesfactor[1,"tau=0"])
}
# NOTE: the 'bayesmeta()' arguments above should probably match
#       the specifications from the original analysis

p <- pppvalue(bma, parameter="tau", statistic=BF, value=0, n=100,
              alternative="greater", rejection.region="lower.tail",
              sigma=bma$sigma)
print(p)

######################################################################
#  If one is only interested in generating samples (and not in test
#  statistics or p-values), one may specify the 'statistic' argument
#  as 'NA'.
#  Note that different 'parameter', 'value' and 'alternative' settings
#  imply different sampling schemes.

p <- pppvalue(bma, parameter="mu", statistic=NA, value=0,
              alternative="less", n=100)

plot(bma, which=2)
abline(h=p$null.value)
points(p$replicates$tau, p$replicates$mu, col="cyan")
}
}
\keyword{ htest }
