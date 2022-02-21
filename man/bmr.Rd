\name{bmr}
\alias{bmr}
\alias{bmr.default}
\alias{bmr.escalc}
\alias{print.bmr}
\alias{plot.bmr}
\alias{pairs.bmr}
\concept{random effects meta analysis}
\concept{random effects meta regression}
\concept{random effects model}
\title{
  Bayesian random-effects meta-regression
}
\description{
  This function allows to derive the posterior distribution of the
  parameters in a random-effects meta-regression and provides functions to
  evaluate joint and marginal posterior probability distributions, etc.
}
\usage{
  bmr(y, ...)
  \method{bmr}{default}(y, sigma, labels = names(y),
    X = matrix(1.0, nrow=length(y), ncol=1,
               dimnames=list(labels,"intercept")),
    tau.prior = "uniform",
    beta.prior.mean = NULL,
    beta.prior.sd = NULL,
    beta.prior.cov = diag(beta.prior.sd^2,
                          nrow=length(beta.prior.sd),
                          ncol=length(beta.prior.sd)),
    interval.type = c("shortest", "central"),
    delta = 0.01, epsilon = 0.0001,
    rel.tol.integrate = 2^16*.Machine$double.eps,
    abs.tol.integrate = 0.0,
    tol.uniroot = rel.tol.integrate, ...)
  \method{bmr}{escalc}(y, labels = NULL, ...)
}
\arguments{
  \item{y}{
    vector of estimates, \emph{or} an \code{\link[metafor]{escalc}} object.
}
  \item{sigma}{
    vector of standard errors associated with \code{y}.
}
  \item{labels}{
    (optional) a vector of labels corresponding to \code{y} and \code{sigma}.
}
  \item{X}{
    (optional) the \emph{regressor matrix} for the regression.
}
  \item{tau.prior}{
    a \code{function} returning the prior density for the heterogeneity
    parameter (\eqn{\tau}) \emph{or} a \code{character} string
    specifying one of the \emph{default \sQuote{non-informative}
    priors}; possible choices for the latter case are: 
    \itemize{
      \item \code{"uniform"}: a uniform prior in \eqn{\tau}
      \item \code{"sqrt"}: a uniform prior in \eqn{\sqrt{\tau}}{sqrt(tau)}
      \item \code{"Jeffreys"}: the Jeffreys prior for \eqn{\tau}
      \item \code{"BergerDeely"}: the prior due to Berger and Deely (1988)
      \item \code{"conventional"}: the conventional prior
      \item \code{"DuMouchel"}: the DuMouchel prior
      \item \code{"shrinkage"}: the \sQuote{uniform shrinkage} prior
      \item \code{"I2"}: a uniform prior on the \sQuote{relative heterogeneity} \eqn{I^2}
    }
    The default is \code{"uniform"} (which should be used with
    caution). The above priors are described in some more detail in the
    \code{\link{bayesmeta}()} help.
}
  \item{beta.prior.mean, beta.prior.sd, beta.prior.cov}{
    the mean and standard deviations, or covariance of the normal prior
    distribution for the effects \eqn{\beta}. If unspecified, an
    (improper) uniform prior is used.
}
  \item{interval.type}{
    the type of (credible, prediction, shrinkage) interval to be
    returned by default; either \code{"shortest"} for shortest
    intervals, or \code{"central"} for central, equal-tailed intervals.
}
  \item{delta, epsilon}{
    the parameters specifying the desired accuracy for approximation of
    the \eqn{\beta} posterior(s), and with that determining the number of
    \eqn{\tau} support points being used internally. Smaller values
    imply greater accuracy and greater computational burden (Roever and
    Friede, 2017).
}
  \item{rel.tol.integrate, abs.tol.integrate, tol.uniroot}{
    the \code{rel.tol}, \code{abs.tol} and \code{tol}
    \sQuote{accuracy} arguments that are passed to
    the \code{\link{integrate}()} or \code{\link{uniroot}()} functions
    for internal numerical integration or root finding
    (see also the help there).
}
  \item{...}{
    other \code{bmr} arguments.
  }
}
\details{
  The random-effects meta-regression model may be stated as
  \deqn{y_i|x_i,\beta,\sigma_i,\tau \;\sim\; \mathrm{Normal}(\beta_1 x_{i,1}
    + \beta_2 x_{i,2} + \ldots + \beta_d x_{i,d}, \;
    \sigma_i^2 + \tau^2)}{y[i] | x[i], beta, sigma[i], tau  ~
    Normal(beta[1]x[i,1] + beta[2]x[i,2] + ... + beta[d]x[i,d], sigma[i]^2 + tau^2)}
  where the data (\eqn{y}, \eqn{\sigma}) enter as \eqn{y_i}{y[i]}, the
  \eqn{i}-th estimate, that is associated with standard error
  \eqn{\sigma_i > 0}{sigma[i] > 0}, where \eqn{i=1,...,k}. In addition to
  estimates and standard errors for the \eqn{i}th observation,
  a set of covariables \eqn{x_{i,j}}{x[i,j]} with \eqn{j=1,...,d} are available
  for each estimate \eqn{y_i}{y[i]}.

  The model includes \eqn{d+1} unknown parameters,
  namely, the \eqn{d} coefficients (\eqn{\beta_1,...,\beta_d}{beta[1],...,beta[d]}), and the heterogeneity
  \eqn{\tau}. Alternatively, the model may also be formulated via an
  intermediate step as
  \deqn{y_i|\theta_i,\sigma_i \;\sim\; \mathrm{Normal}(\theta_i, \, \sigma_i^2),}{y[i]|theta[i],sigma[i]  ~  Normal(theta[i], sigma[i]^2),}
  \deqn{\theta_i|\beta,x_i,\tau \;\sim\; \mathrm{Normal}(\beta_1 x_{i,1}
  + \ldots + \beta_d x_{i,d}, \; \tau^2),}{theta[i]|beta,x[i],tau  ~  Normal(beta[1]x[i,1]+...+beta[d]x[i,d], tau^2),}
  where the \eqn{\theta_i}{theta[i]} denote the \emph{trial-specific} means
  that are then measured through the estimate \eqn{y_i}{y[i]} with an
  associated measurement uncertainty \eqn{\sigma_i}{sigma[i]}. The
  \eqn{\theta_i}{theta[i]} again differ from trial to trial (even for
  identical covariable vectors \eqn{x_i}{x[i]}) and are
  distributed around a mean of \eqn{\beta_1 x_{i,1}
  + \ldots + \beta_d x_{i,d}}{beta[1]x[i,1]+...+beta[d]x[i,d]} with
  standard deviation \eqn{\tau}.

  It if often convenient to express the model in matrix notation, i.e.,
  \deqn{y|\theta,\sigma \;\sim\; \mathrm{Normal}(\theta, \,
    \Sigma)}{y | theta, sigma  ~
    Normal(theta, Sigma)}
  \deqn{\theta|X,\beta,\tau \;\sim\; \mathrm{Normal}(X \beta, \,
    \tau I)}{theta | X, beta, tau  ~
    Normal(X beta, tau I)}
  where \eqn{y}, \eqn{\sigma}, \eqn{\beta} and \eqn{\theta} now denote
  \eqn{k}-dimensional vectors, \eqn{X} is the (\eqn{k \times d}{k x
  d}) \emph{regressor matrix}, and \eqn{\Sigma} is a (\eqn{k \times k}{k x
  k}) diagonal covariance matrix containing the
  \eqn{\sigma_i^2}{sigma[i]^2} values, while
  \eqn{I} is the (\eqn{k \times k}{k x k}) identity matrix. The
  regressor matrix \eqn{X} plays a crucial role here, as the
  \sQuote{\code{X}} argument (with rows corresponding to studies, and
  columns corresponding to covariables) is required to specify the exact
  regression setup.

  Meta-regression allows the consideration of (study-level) covariables
  in a meta-analysis. Quite often, these may also be indicator variables
  (\sQuote{zero/one} variables) simply identifying subgroups of studies.
  See also the examples shown below.

  \subsection{Connection to the simple random-effects model}{
  The meta-regression model is a generalisation of the \sQuote{simple}
  random-effects model that is implemented in the
  \code{\link{bayesmeta}()} function. Meta-regression reduces to the
  estimation of a single \dQuote{intercept} term when the regressor
  matrix (\eqn{X}) consists of a single column of
  ones (which is also the default setting in case the \sQuote{\code{X}}
  argument is left unspecified). The single regression coefficient
  \eqn{\beta_1}{beta[1]} then is equivalent to the \eqn{\mu} parameter
  from the simple random effects model (see also the \sQuote{Examples}
  section below).
  }

  \subsection{Specification of the regressor matrix}{ The actual
   regression model is specified through the regressor matrix \eqn{X},
   which is supplied via the \sQuote{\code{X}} argument, and which
   often may be specified in different ways. There usually is no unique
   solution, and what serves the present purpose best then depends on
   the context; see also the examples below. Sensible column names
   should be specified for \code{X}, as these will subsequently
   determine the labels for the associated parameters later on. Model
   specification via the regressor matrix has the advantage of being
   very explicit and transparent; if one prefers a
   code{\link[stats]{formula}} interface instead, a regressor matrix may
   be generated via the \sQuote{\code{\link[stats]{model.matrix}()}}
   function.
  }

  \subsection{Prior specification}{
  Priors for \eqn{\beta} and \eqn{\tau} are assumed to factor into
  into independent marginals \eqn{p(\beta,\tau)=p(\beta)\times
  p(\tau)}{p(beta, tau) = p(beta) * p(tau)} and either (improper)
  uniform or a normal priors may be specified for the regression coefficients
  \eqn{\beta}. 
  For sensible prior choices for the heterogeneity parameter \eqn{\tau},
  see also Roever (2020), Roever \emph{et al.} (2021) and the
  \sQuote{\code{\link{bayesmeta}()}} function's help.
  }
  
  \subsection{Accessing posterior density functions, etc.}{
    Within the \code{\link{bayesmeta}()} function, access to posterior
    density, cumulative distribution function, quantile functtion,
    random number generation and posterior inverval computation is
    implemented via the \code{$dposterior()}, \code{$dposterior()},
    \code{$pposterior()}, \code{$qposterior()}, \code{$rposterior()}
    and \code{$post.interval()} functions that are accessible as elements
    in the returned \code{list} object. Prediction and shrinkage
    estimation are available by setting additional arguments in the
    above functions.
    
    In the meta-regression context things get slightly more complicated,
    as the \eqn{\beta} parameter may be of higher dimension. Hence, in the
    \code{bmr()} function, the three different types of distributions
    related to \emph{posterior distribution}, \emph{prediction} and
    \emph{shrinkage} are split up into three groups of
    functions. For example, the posterior density is accessible via the
    \code{$dposterior()} function, the predictive distribution via the
    \code{$dpredict()} function, and the shrinkage estimates via the
    \code{$dshrink()} function. Analogous functions are returned for
    cumulative distribution, quantile function, etc.; see also the
    \sQuote{Value} section below.
  }
    
  \subsection{Computation}{
  The \code{bmr()} function utilizes the same computational method
  as the \code{\link{bayesmeta}()} function to derive the posterior
  distribution, namely, the \acronym{DIRECT} algorithm. Numerical
  accuracy of the computations is determined by the \sQuote{\code{delta}}
  and \sQuote{\code{epsilon}} arguments (Roever and Friede,
  2017).

  A slight difference between the \code{\link{bayesmeta}()} and
  \code{bmr()} implementations exists in the determination of the grid
  approximation within the \acronym{DIRECT} algorithm. While
  \code{bmr()} considers divergences w.r.t. the conditional posterior
  distributions \eqn{p(\beta|\tau)}, \code{bayesmeta()} in addition
  considers divergences w.r.t. the shrinkage estimates, which in general
  leads to a denser binning (as one can see from the numbers of bins
  required; see the example below). A denser binning within the
  \code{bmr()} function may be achieved by reducing the
  \sQuote{\code{delta}} argument.
  }
}
\value{
  A list of class \code{bmr} containing the following elements:
  \item{y}{vector of estimates (the input data).}
  \item{sigma}{vector of standard errors corresponding
    to \code{y} (input data).}
  \item{X}{the regressor matrix.}
  \item{k}{number of data points (length of \code{y}, or rows of \eqn{X}).}
  \item{d}{number of coefficients (columns of \code{X}).}
  \item{labels}{vector of labels corresponding to \code{y} and \code{sigma}.}
  \item{variables}{variable names for the \eqn{\beta} coefficients
  (determined by the column names of the supplied \code{X} argument).}
  \item{tau.prior}{the prior probability density function for \eqn{\tau}.}
  \item{tau.prior.proper}{a \code{logical} flag indicating whether the
    heterogeneity prior appears to be proper (which is judged based on
    an attempted numerical integration of the density function).} 
  \item{beta.prior}{a \code{list} containing the prior mean vector and
  covariance matrix for the coefficients \eqn{\beta}.}
  \item{beta.prior.proper}{a \code{logical} vector (of length \eqn{d})
  indicating whether the corresponding \eqn{\beta} coefficient's prior is
  proper (i.e., finite prior mean and variance were specified).} 
  \item{dprior}{a \code{function(tau, beta, which.beta, log=FALSE)} of
  \eqn{\tau} and/or \eqn{\beta} parameters, returning either the joint or
    marginal prior probability density, depending on which parameter(s)
    is/are provided.}
  \item{likelihood}{a \code{function(tau, beta, which.beta)} \eqn{\tau}
  and/or \eqn{\beta}, returning either the joint or marginal
  likelihood, depending on which parameter(s) is/are provided.}
  \item{dposterior, pposterior, qposterior, rposterior, post.interval}{functions of \eqn{\tau} and/or \eqn{\beta} parameters,
    returning either the joint or marginal posterior probability density,
    (depending on which parameter(s) is/are provided), or cumulative
  distribution function, quantile function, random numbers or posterior intervals.}
  \item{dpredict, ppredict, qpredict, rpredict, pred.interval}{functions
  of \eqn{\beta} returning density, cumulative distribution
    function, quantiles, random numbers, or intervals for the
  \emph{predictive distribution}. This requires specification of \eqn{x}
  values to indicate what covariable values to consider. Use of
  \sQuote{\code{mean=TRUE}} (the default) yields predictions for the
  \emph{mean} (\eqn{x'\beta}{t(x) beta} values), setting it to
  \code{FALSE} yields \emph{predictions} (\eqn{\theta} values).}
  \item{dshrink, pshrink, qshrink, rshrink, shrink.interval}{functions
  of \eqn{\theta} yielding density, cumulative distribution, quantiles,
  random numbers or posterior intervals for the \emph{shrinkage
  estimates} of the individual \eqn{\theta_i}{theta[i]} parameters
  corresponding to the supplied \eqn{y_i}{y[i]} data values
  (\eqn{i=1,\ldots,k}{i=1,...,k}). May be identified using the
  \sQuote{\code{which}} argument via its index (\eqn{i}) or a character
  string giving the corresponding study label.}
  \item{post.moments}{a \code{function(tau)}
    returning conditional posterior moments (mean and covariance) of
    \eqn{\beta} as a function of \eqn{\tau}.}
  \item{pred.moments}{a \code{function(tau, x, mean=TRUE)}
    returning conditional posterior predictive moments (means and
    standard deviations) as a function of \eqn{\tau}.}
  \item{shrink.moments}{a \code{function(tau, which)}
    returning conditional moments (means and standard deviations of
    shrinkage distributions) as a function of \eqn{\tau}.}
  \item{summary}{a \code{matrix} listing some summary statistics, namely
    marginal posterior mode, median, mean, standard deviation
    and a (shortest) 95\% credible intervals,
    of the marginal posterior distributions of \eqn{\tau} and \eqn{\beta_i}{beta[i]}.}
  \item{interval.type}{the \code{interval.type} input argument
    specifying the type of interval to be returned by default.}
  \item{ML}{a \code{matrix} giving joint and marginal maximum-likelihood
    estimates of \eqn{(\tau,\beta)}.}
  \item{MAP}{a \code{matrix} giving joint and marginal
    maximum-a-posteriori estimates of \eqn{(\tau,\beta)}.}
  \item{theta}{a \code{matrix} giving the \sQuote{shrinkage estimates},
    i.e, summary statistics of the trial-specific means
    \eqn{\theta_i}{theta[i]}.} 
  \item{marginal.likelihood}{the marginal likelihood of the data (this 
    number can only be computed if proper effect and heterogeneity priors 
    are specified).}
  \item{bayesfactor}{Bayes factors and minimum Bayes factors for the
    hypotheses of \eqn{\tau=0} and \eqn{\beta_i=0}{beta[i]=0}. These depend on the
    marginal likelihood and hence can only be computed if proper effect
    and/or heterogeneity priors are specified.}
  \item{support}{a \code{list} giving the \eqn{\tau} support points used
    internally in the grid approximation, along with their associated
    weights, and conditional mean and covariance of \eqn{\beta}.}
  \item{delta, epsilon}{the \sQuote{\code{delta}} and \sQuote{\code{epsilon}}
    input parameter determining numerical accuracy.}
  \item{rel.tol.integrate, abs.tol.integrate, tol.uniroot}{the input
    parameters determining the numerical accuracy of the internally used
    \code{\link{integrate}()} and \code{\link{uniroot}()} functions.}
  \item{call}{an object of class \code{call} giving the
    function call that generated the \code{bmr} object.}
  \item{init.time}{the computation time (in seconds) used to generate
    the \code{bmr} object.}
  }
\references{
  C. Roever.
  Bayesian random-effects meta-analysis using the bayesmeta R package.
  \emph{Journal of Statistical Software}, \bold{93}(6):1-51, 2020.
  \doi{10.18637/jss.v093.i06}.

  C. Roever, R. Bender, S. Dias, C.H. Schmid, H. Schmidli, S. Sturtz,
  S. Weber, T. Friede.
  On weakly informative prior distributions for the heterogeneity
  parameter in Bayesian random-effects meta-analysis.
  \emph{Research Synthesis Methods}, \bold{12}(4):448-474, 2021.
  \doi{10.1002/jrsm.1475}.
  
  C. Roever, T. Friede.
  Discrete approximation of a mixture distribution via restricted divergence.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{26}(1):217-222, 2017.
  \doi{10.1080/10618600.2016.1276840}.

  A. Gelman, J.B. Carlin, H.S. Stern, D.B. Rubin. \emph{Bayesian data
  analysis}. Chapman & Hall / CRC, Boca Raton, 1997.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{bayesmeta}}, \code{\link[metafor]{escalc}},
  \code{\link[stats]{model.matrix}}, \code{\link{CrinsEtAl2014}},
  \code{\link{RobergeEtAl2017}}.
}
\examples{
\dontrun{
######################################################################
# (1)  A simple example with two groups of studies

# load data:
data("CrinsEtAl2014")
# compute effect measures (log-OR):
crins.es <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)
# show data:
crins.es[,c("publication", "IL2RA", "exp.AR.events", "exp.total",
            "cont.AR.events", "cont.total", "yi", "vi")]

# specify regressor matrix:
X <- cbind("bas"=as.numeric(crins.es$IL2RA=="basiliximab"),
           "dac"=as.numeric(crins.es$IL2RA=="daclizumab"))
print(X)
print(cbind(crins.es[,c("publication", "IL2RA")], X))
# NB: regressor matrix specifies individual indicator covariates
#     for studies with "basiliximab" and "daclizumab" treatment.

# perform regression:
bmr01 <- bmr(y=crins.es$yi, sigma=sqrt(crins.es$vi),
             labels=crins.es$publication, X=X)

# alternatively, one may simply supply the "escalc" object
# (yields identical results):
bmr01 <- bmr(crins.es, X=X)

# show results:
bmr01
bmr01$summary
plot(bmr01)
pairs(bmr01)

# NOTE: there are many ways to set up the regressor matrix "X"
# (also affecting the interpretation of the involved parameters).
# See the above specification and check out the following alternatives:
X <- cbind("bas"=1, "offset.dac"=c(1,0,1,0,0,0))
X <- cbind("intercept"=1, "offset"=0.5*c(1,-1,1,-1,-1,-1))
# One may also use the "model.matrix()" function
# to specify regressor matrices via the "formula" interface; e.g.:
X <- model.matrix( ~ IL2RA, data=crins.es)
X <- model.matrix( ~ IL2RA - 1, data=crins.es)


######################################################################
# (2)  A simple example reproducing a "bayesmeta" analysis:

data("CrinsEtAl2014")
crins.es <- escalc(measure="OR",
                   ai=exp.AR.events,  n1i=exp.total,
                   ci=cont.AR.events, n2i=cont.total,
                   slab=publication, data=CrinsEtAl2014)

# a "simple" meta-analysis:
bma02 <- bayesmeta(crins.es,
                   tau.prior=function(t){dhalfnormal(t, scale=0.5)},
                   mu.prior.mean=0, mu.prior.sd=4)

# the equivalent "intercept-only" meta-regression:
bmr02 <- bmr(crins.es,
             tau.prior=function(t){dhalfnormal(t, scale=0.5)},
             beta.prior.mean=0, beta.prior.sd=4)
# the corresponding (default) regressor matrix:
bmr02$X

# compare computation time and numbers of bins used internally:
cbind("seconds" = c("bayesmeta" = unname(bma02$init.time),
                    "bmr"       = unname(bmr02$init.time)),
      "bins"    = c("bayesmeta" = nrow(bma02$support),
                    "bmr"       = nrow(bmr02$support$tau)))

# compare heterogeneity estimates:
rbind("bayesmeta"=bma02$summary[,1], "bmr"=bmr02$summary[,1])

# compare effect estimates:
rbind("bayesmeta"=bma02$summary[,2], "bmr"=bmr02$summary[,2])


######################################################################
# (3)  An example with binary as well as continuous covariables:

# load data:
data("RobergeEtAl2017")
str(RobergeEtAl2017)
head(RobergeEtAl2017)
?RobergeEtAl2017

# compute effect sizes (log odds ratios) from count data:
es.pe  <- escalc(measure="OR",
                 ai=asp.PE.events,  n1i=asp.PE.total,
                 ci=cont.PE.events, n2i=cont.PE.total,
                 slab=study, data=RobergeEtAl2017,
                 subset=complete.cases(RobergeEtAl2017[,7:10]))

# show "bubble plot" (bubble sizes are
# inversely proportional to standard errors):
plot(es.pe$dose, es.pe$yi, cex=1/sqrt(es.pe$vi),
     col=c("blue","red")[as.numeric(es.pe$onset)],
     xlab="dose (mg)", ylab="log-OR (PE)", main="Roberge et al. (2017)")
legend("topright", col=c("blue","red"), c("early onset", "late onset"), pch=1)

# set up regressor matrix:
# (individual intercepts and slopes for two subgroups):
X <- model.matrix(~ -1 + onset + onset:dose, data=es.pe)
colnames(X) <- c("intEarly", "intLate", "slopeEarly", "slopeLate")
# check out regressor matrix (and compare to original data):
print(X)

# perform regression:
bmr03 <- bmr(es.pe, X=X)
bmr03$summary

# derive predictions from the model;
# specify corresponding "regressor matrices":
newx.early <- cbind(1, 0, seq(50, 150, by=5), 0)
newx.late  <- cbind(0, 1, 0, seq(50, 150, by=5))
# (note: columns correspond to "beta" parameters)

# compute predicted medians and 95 percent intervals: 
pred.early <- cbind("median"=bmr03$qpred(0.5, x=newx.early),
                    bmr03$pred.interval(x=newx.early))
pred.late <- cbind("median"=bmr03$qpred(0.5, x=newx.late),
                    bmr03$pred.interval(x=newx.late))

# draw "bubble plot": 
plot(es.pe$dose, es.pe$yi, cex=1/sqrt(es.pe$vi),
     col=c("blue","red")[as.numeric(es.pe$onset)],
     xlab="dose (mg)", ylab="log-OR (PE)", main="Roberge et al. (2017)")
legend("topright", col=c("blue","red"), c("early onset", "late onset"), pch=1)
# add predictions to bubble plot:
matlines(newx.early[,3], pred.early, col="blue", lty=c(1,2,2))
matlines(newx.late[,4], pred.late, col="red", lty=c(1,2,2))

}
}
\keyword{ models }
