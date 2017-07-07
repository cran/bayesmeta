\name{bayesmeta}
\alias{bayesmeta}
\alias{bayesmeta.default}
\alias{bayesmeta.escalc}
\alias{print.bayesmeta}
\alias{summary.bayesmeta}
\concept{random effects meta analysis}
\concept{random effects model}
\title{
  Bayesian random-effects meta-analysis
}
\description{
  This function allows to derive the posterior distribution of the two
  parameters in a random-effects meta-analysis and provides functions to
  evaluate joint and marginal posterior probability distributions, etc.
}
\usage{
  bayesmeta(y, ...)
  \method{bayesmeta}{default}(y, sigma, labels = names(y),
          tau.prior = "uniform",
          mu.prior = c("mean"=NA,"sd"=NA),
          mu.prior.mean = mu.prior[1], mu.prior.sd = mu.prior[2],
          interval.type = c("shortest", "central"),
          delta = 0.01, epsilon = 0.0001,
          rel.tol.integrate = 2^16*.Machine$double.eps,
          abs.tol.integrate = rel.tol.integrate, ...)
  \method{bayesmeta}{escalc}(y, labels = NULL, ...)
}
\arguments{
  \item{y}{
    vector of estimates \emph{or} an \code{\link[metafor]{escalc}} object.
}
  \item{sigma}{
    vector of standard errors associated with \code{y}.
}
  \item{labels}{
    (optional) a vector of labels corresponding to \code{y} and \code{sigma}.
}
  \item{tau.prior}{
    a \code{function} returning the prior density for the heterogeneity
    parameter (\eqn{\tau}) \emph{or} a \code{character} string
    specifying one of the \emph{default \sQuote{non-informative}
    priors}; possible choices for the latter case are: 
    \itemize{
      \item \code{"uniform"}: a uniform prior in \eqn{\tau}
      \item \code{"sqrt"}: a uniform prior in \eqn{\sqrt{\tau}}{sqrt(tau)}.
      \item \code{"Jeffreys"}: the Jeffreys prior for \eqn{\tau}
      \item \code{"BergerDeely"}: the prior due to Berger and Deely (1988)
      \item \code{"DuMouchel"}: the DuMouchel prior
      \item \code{"shrinkage"}: the \sQuote{uniform shrinkage} prior
      \item \code{"I2"}: a uniform prior on the \sQuote{relative heterogeneity} \eqn{I^2}
    }
    The default is \code{"uniform"} (which should be used with caution;
    see remarks below).
}
  \item{mu.prior}{
    the mean and standard deviation of the normal prior distribution for
    the effect \eqn{\mu}. If unspecified, an (improper) uniform prior is
    used.
}
  \item{mu.prior.mean, mu.prior.sd}{
    alternative parameters to specify the prior distribution for the
    effect \eqn{\mu} (see above).
}
  \item{interval.type}{
    the type of (credibility, prediction, shrinkage) interval to be
    returned by default; either \code{"shortest"} for shortest
    intervals, or \code{"central"} for central, equal-tailed intervals.
}
  \item{delta, epsilon}{
    the parameters specifying the desired accuracy for approximation of
    the \eqn{\mu} posterior(s), and with that determining the number of
    \eqn{\tau} support points being used internally. Smaller values
    imply greater accuracy and greater computational burden (Roever and
    Friede, 2017).
}
  \item{rel.tol.integrate, abs.tol.integrate}{
    the \code{rel.tol} and \code{abs.tol} arguments that are passed to
    the \code{\link{integrate}()} function for internal numerical integration
    (see also the help there).
}
  \item{...}{
    other \code{bayesmeta} arguments.
}
}
\details{
  The common random-effects meta-analysis model may be stated as
  \deqn{y_i|\mu,\sigma_i,\tau \;\sim\; \mathrm{Normal}(\mu, \, \sigma_i^2 + \tau^2)}{y[i] | mu, sigma[i], tau  ~  Normal(mu, sigma[i]^2 + tau^2)}
  where the data (\eqn{y}, \eqn{\sigma}) enter as \eqn{y_i}{y[i]}, the
  \eqn{i}-th estimate, that is associated with standard error
  \eqn{\sigma_i}{sigma[i]}, and \eqn{i=1,...,k}. The model includes two unknown parameters,
  namely the (mean) effect \eqn{\mu}, and the heterogeneity
  \eqn{\tau}. Alternatively, the model may also be formulated via an
  intermediate step as
  \deqn{y_i|\theta_i,\sigma_i \;\sim\; \mathrm{Normal}(\theta_i, \, \sigma_i^2),}{y[i]|theta[i],sigma[i]  ~  Normal(theta[i], sigma[i]^2),}
  \deqn{\theta_i|\mu,\tau \;\sim\; \mathrm{Normal}(\mu, \, \tau^2),}{theta[i]|mu,tau  ~  Normal(mu, tau^2),}
  where the \eqn{\theta_i}{theta[i]} denote the \emph{trial-specific} means
  that are then measured through the estimate \eqn{y_i}{y[i]} with an
  associated measurement uncertainty \eqn{\sigma_i}{sigma[i]}. The
  \eqn{\theta_i}{theta[i]} again differ from trial to trial and are
  distributed around a common mean \eqn{\mu} with standard deviation
  \eqn{\tau}.

  The \code{bayesmeta()} function utilizes the fact that the joint posterior
  distribution \eqn{p(\mu, \tau | y, \sigma)} may be partly analytically
  integrated out to determine the integrals necessary for coherent
  Bayesian inference on one or both of the parameters.
  
  As long as we assume that the prior probability distribution may be
  factored into independent marginals \eqn{p(\mu,\tau)=p(\mu)\times
  p(\tau)}{p(mu, tau) = p(mu) * p(tau)} and either an (improper) uniform
  or a normal prior is used for the effect \eqn{\mu}, the joint
  likelihood \eqn{p(y|\mu,\tau)} may be analytically marginalized over
  \eqn{\mu}, yielding the marginal likelihood function
  \eqn{p(y|\tau)}. Using any prior \eqn{p(\tau)} for the heterogeneity,
  the 1-dimensional marginal posterior \eqn{p(\tau|y) \propto
  p(y|\tau)\times p(\tau)}{p(tau|y) = p(y|tau) * p(tau) * const.} may
  then be treated numerically. As the \emph{conditional} posterior
  \eqn{p(\mu|\tau,y)} is a normal distribution, inference on the
  remaining joint (\eqn{p(\mu,\tau|y)}) or marginal (\eqn{p(\mu|y)})
  posterior may be approached numerically as well. Accuracy of the
  computations is determined by the \code{delta} (maximum divergence
  \eqn{\delta}) and \code{epsilon} (tail probability \eqn{\epsilon})
  parameters (Roever and Friede, 2017).

  What constitutes a sensible prior for both \eqn{\mu} and \eqn{\tau}
  depends (as usual) very much on the context.
  Potential candidates for informative (or weakly informative)
  heterogeneity (\eqn{\tau}) priors may include the \emph{half-normal},
  \emph{half-Student-\eqn{t}},  \emph{half-Cauchy}, \emph{exponential},
  or \emph{Lomax} distributions; we recommend the use of heavy-tailed
  prior distributions if in case of prior/data conflict the prior is
  supposed to be discounted (O'Hagan and Pericchi, 2012).  A sensible
  informative prior might also be a \emph{log-normal} or a \emph{scaled
  inverse \eqn{\chi^2}} distribution. 
  One might argue that an uninformative prior for \eqn{\tau} may be
  uniform or monotonically decreasing in \eqn{\tau}. Another option is
  to use the \emph{Jeffreys prior} (see the \code{tau.prior} argument
  above). The Jeffreys prior implemented here is the variant also
  described by Tibshirani (1989) that results from fixing the location
  parameter (\eqn{\mu}) and considering the Fisher information element
  corresponding to the heterogeneity \eqn{\tau} only. This prior also
  constitutes the \emph{Berger/Bernardo reference prior} for the present
  problem (Bodnar et al., 2016).  The \emph{uniform shrinkage} and
  \emph{DuMouchel} priors are described in Spiegelhalter et al. (2004,
  Sec. 5.7.3). 
  The procedure is able to handle improper priors (and does so by
  default), but of course the usual care must be taken here, as the
  resulting posterior \emph{may} or \emph{may not} be a proper
  probability distribution.
  
  Note that a wide range of different types of endpoints may be treated
  on a continuous scale after an appropriate transformation; for
  example, count data may be handled by considering corresponding
  logarithmic odds ratios. Many such transformations are implemented
  in the \pkg{metafor} package's \code{\link[metafor]{escalc}} function
  and may be directly used as an input to the \code{bayesmeta()}
  function; see also the example below.  Alternatively, the
  \pkg{compute.es} package also provides a range of effect sizes to be
  computed from different data types.

  The \code{bayesmeta()} function eventually generates some basic
  summary statistics, but most importantly it provides an object
  containing a range of \code{function}s allowing to evaluate posterior
  distributions; this includes joint and marginal posteriors, prior and
  likelihood, predictive distributions, densities, cumulative
  distributions and quantile functions. For more details see also the
  documentation and examples below.
  Use of the \code{individual} argument allows to access posteriors
  of study-specific (\emph{shrinkage-}) effects
  (\eqn{\theta_i}{theta[i]}).
  The \code{predict} argument may be used to access the predictive
  distribution of a future study's effect
  (\eqn{\theta_{k+1}}{theta[k+1]}), facilitating a
  \emph{meta-analytic-predictive (MAP)} approach (Neuenschwander et al.,
  2010).

  \subsection{Prior specification details}{
    When specifying the \code{tau.prior} argument as a \code{character} string
    (and not as a prior density \code{function}), then the actual
    prior probability density functions corresponding to the possible
    choices of the \code{tau.prior} argument are given by:
    \itemize{
      \item \code{"uniform"} - the (improper) uniform prior in \eqn{\tau}:
        \deqn{p(\tau) \;\propto\; 1}{p(tau) = 1}
      \item \code{"sqrt"} - the (improper) uniform prior in \eqn{\sqrt{\tau}}{sqrt(tau)}:
        \deqn{p(\tau) \;\propto\; \tau^{-1/2} \;=\; \frac{1}{\sqrt{\tau}}}{p(tau) = tau^(-1/2) = 1/sqrt(tau)}
      \item \code{"Jeffreys"} - \emph{Tibshirani's noninformative prior},
        a variation of the \emph{Jeffreys prior}, which here also constitutes
        the \emph{Berger/Bernardo reference prior} for \eqn{\tau}:
        \deqn{p(\tau) \;\propto\;
        \sqrt{\sum_{i=1}^k\Bigl(\frac{\tau}{\sigma_i^2+\tau^2}\Bigr)^2}}{p(tau) = sqrt(sum((tau/(sigma[i]^2+tau^2))^2))}
        (This is also an improper prior whose density does not integrate
  to 1).
      \item \code{"BergerDeely"} - the (improper) \emph{Berger/Deely} prior:
        \deqn{p(\tau) \;\propto\; \prod_{i=1}^k \Bigl(\frac{\tau}{\sigma_i^2+\tau^2}\Bigr)^{1/k}}{p(tau) = prod((tau / (sigma[i]^2+tau^2))^(1/k))}
        This is a variation of the above \emph{Jeffreys} prior, and both are equal in
        case all standard errors (\eqn{\sigma_i}{sigma[i]}) are the same.
      \item \code{"DuMouchel"} - the (proper) \emph{DuMouchel} prior:
        \deqn{p(\tau) \;=\; \frac{s_0}{(s_0+\tau)^2}}{p(tau) = s0 / (s0+tau)^2}
        where \eqn{s_0=\sqrt{s_0^2}}{s0=sqrt(s0^2)} and \eqn{s_0^2}{s0^2} again is the harmonic mean of the standard errors (as above).
      \item \code{"shrinkage"} - the (proper) \emph{uniform shrinkage} prior:
        \deqn{p(\tau) \;=\; \frac{2 s_0^2 \tau}{(s_0^2+\tau^2)^2}}{p(tau)=(2*s0^2*tau) / (s0^2+tau^2)^2}
        where \eqn{s_0^2=\frac{k}{\sum_{i=1}^k \sigma_i^{-2}}}{s0^2 = k/sum(sigma[i]^(-2))} is the harmonic
        mean of the squared standard errors \eqn{\sigma_i^2}{sigma[i]^2}.
      \item \code{"I2"} - the (proper) uniform prior in \eqn{I^2}:
        \deqn{p(\tau) \;=\; \frac{2 \hat{\sigma}^2 \tau}{(\hat{\sigma}^2 + \tau^2)^2}}{p(tau)=(2*sigmaHat^2*tau) /(sigmaHat^2+tau^2)^2}
        where \eqn{\hat{\sigma}^2 = \frac{(k-1)\; \sum_{i=1}^k\sigma_i^{-2}}{\bigl(\sum_{i=1}^k\sigma_i^{-2}\bigr)^2 - \sum_{i=1}^k\sigma_i^{-4}}}{sigmaHat^2 = (k-1)*sum(1/sigma[i]^2) / (sum(1/sigma[i]^2)^2 - sum(1/sigma[i]^4))}.
        This prior is similar to the uniform shrinkage prior, except for
        the use of \eqn{\hat{\sigma}^2}{sigmaHat^2} instead of \eqn{s_0^2}{s0^2}.
    }
  }

  \subsection{Credibility intervals}{
    Credibility intervals (as well as prediction and shrinkage intervals)
    may be determined in different ways. By default, \emph{shortest}
    intervals are returned, which for unimodal posteriors (the usual
    case) is equivalent to the \emph{highest posterior density region}
    (Gelman et al., 1997, Sec. 2.3).
    Alternatively, central (equal-tailed) intervals may also be derived.
    The default behaviour may be controlled via the \code{interval.type}
    argument, or also by using the \code{method} argument with each
    individual call of the \code{$post.interval()} function (see below).
    A third option, although not available for prediction or shrinkage
    intervals, and hence not as an overall default choice, but only for
    the \code{$post.interval()} function, is to
    determine the \emph{evidentiary} credibility interval, which has the
    advantage of being parameterization invariant (Shalloway, 2014).
  }
  
  \subsection{Bayes factors}{
    Bayes factors (Kass and Raftery, 1995) for the two hypotheses of
    \eqn{\tau=0} and \eqn{\mu=0} are provided in the \code{$bayesfactor}
    element; \emph{low} or \emph{high} values here constitute evidence
    \emph{against} or \emph{in favour of} the hypotheses,
    respectively. Bayes factors are based on marginal likelihoods and
    can only be computed if the priors for heterogeneity and effect are
    proper. Bayes factors for other hypotheses can be computed using the
    marginal likelihood (as provided through the \code{$marginal}
    element) and the \code{$likelihood()} function. Bayes factors must
    be interpreted with very much caution, as they are susceptible to
    \emph{Lindley's paradox} (Lindley, 1957), which especially implies
    that variations of the prior specifications that have only minuscule
    effects on the posterior distribution may have a substantial impact
    on Bayes factors (via the marginal likelihood). For more details on
    the problems and challenges related to Bayes factors see also
    Gelman et al. (1997, Sec. 7.4).
  }
  
  \subsection{Numerical accuracy}{
    Accuracy of the numerical results is determined by four parameters,
    namely, the accuracy of numerical integration as specified through the
    \code{rel.tol.integrate} and \code{abs.tol.integrate} arguments (which
    are internally passed on to the \code{\link[stats]{integrate}}
    function), and the accuracy of the grid approximation used for
    integrating out the heterogeneity as specified through the
    \code{delta} and \code{epsilon} arguments (Roever and Friede,
    2017). As these may also heavily impact on the computation time, be
    careful when changing these from their default values. You can monitor
    the effect of different settings by checking the number and range of
    support points returned in the \code{$support} element.
  }  
}
\value{
  A list of class \code{bayesmeta} containing the following elements:
  \item{y}{vector of estimates (the input data).}
  \item{sigma}{vector of standard errors corresponding
    to \code{y} (input data).}
  \item{labels}{vector of labels corresponding to \code{y} and \code{sigma}.}
  \item{k}{number of data points (in \code{y}).}
  \item{tau.prior}{the prior probability density function for \eqn{\tau}.}
  \item{mu.prior.mean}{the prior mean of \eqn{\mu}.}
  \item{mu.prior.sd}{the prior standard deviation of \eqn{\mu}.}
  \item{dprior}{a \code{function(tau=NA, mu=NA, log=FALSE)} of two
    parameters, \code{tau} and/or \code{mu}, returning either the joint or
    marginal prior probability density, depending on which parameter(s)
    is/are provided.} 
  \item{tau.prior.proper}{a \code{logical} flag indicating whether the
    heterogeneity prior appears to be proper (which is judged based on
    an attempted numerical integration of the density function).} 
  \item{likelihood}{a \code{function(tau=NA, mu=NA, log=FALSE)} of two
    parameters, \code{tau} and/or \code{mu}, returning either the joint or
    marginal likelihood, depending on which parameter(s) is/are provided.} 
  \item{dposterior}{a \code{function(tau=NA, mu=NA, theta=mu, log=FALSE, predict=FALSE, individual=FALSE)}
    of two parameters, \code{tau} and/or \code{mu},
    returning either the joint or marginal posterior probability density,
    depending on which parameter(s) is/are provided. Using the argument
    \code{predict=TRUE} yields the \emph{posterior predictive
    distribution} for \eqn{\theta}. Using the \code{individual}
    argument, you can request individual effects'
    (\eqn{\theta_i}{theta[i]}) posterior distributions. May be an integer
    number (\code{1,...,k}) giving the index, or a character string giving
    the label.}
  \item{pposterior}{a \code{function(tau=NA, mu=NA, theta=mu, predict=FALSE, individual=FALSE)} 
    of one parameter (either \code{tau} or \code{mu}) returning the
    corresponding marginal posterior cumulative distribution
    function. Using the argument \code{predict=TRUE} yields the posterior
    predictive distribution for \eqn{\theta}. Using the \code{individual}
    argument, you can request individual effects'
    (\eqn{\theta_i}{theta[i]}) posterior distributions. May be an integer
    number (\code{1,...,k}) giving the index, or a character string giving
    the label.} 
  \item{qposterior}{a \code{function(tau.p=NA, mu.p=NA, theta.p=mu.p, predict=FALSE, individual=FALSE)}
    of one parameter (either \code{tau.p} or \code{mu.p}) returning the
    corresponding marginal posterior quantile function. Using the argument
    \code{predict=TRUE} yields the posterior predictive distribution for
    \eqn{\theta}. Using the \code{individual} argument, you can request
    individual effects' (\eqn{\theta_i}{theta[i]}) posterior
    distributions. May be an integer number (\code{1,...,k}) giving the
    index, or a character string giving the label.}
  \item{rposterior}{a \code{function(n=1, predict=FALSE, individual=FALSE, tau.sample=TRUE)}
    generating \code{n} independent random draws from the
    (2-dimensional) posterior distribution. Using the argument
    \code{predict=TRUE} yields the posterior predictive distribution for
    \eqn{\theta}. Using the \code{individual} argument, you can request
    individual effects' (\eqn{\theta_i}{theta[i]}) posterior
    distributions. May be an integer number (\code{1,...,k}) giving the
    index, or a character string giving the label. In general, this via
    the inversion method, so it is rather slow. However, if one is not
    interested in sampling the heterogeneity parameter (\eqn{\tau}),
    using \sQuote{\code{tau.sample=FALSE}} will speed up the function
    substantially.} 
  \item{post.interval}{a \code{function(tau.level=NA, mu.level=NA, theta.level=mu.level, method=c("shortest","central","evidentiary"), predict=FALSE, individual=FALSE)}
    returning a credibility interval, depending on which of the two
    parameters is provided (either \code{tau.level} or
    \code{mu.level}). The additional parameter \code{method} may be used
    to specify the desired type of interval: \code{method} \code{=}
    \code{"shortest"} returns the shortest interval, \code{method}
    \code{=} \code{"central"} returns a central interval, and \code{method}
    \code{=} \code{"evidentiary"} returns an evidentiary interval
    (Shalloway, 2014); the former is the default option. Using the
    argument \code{predict=TRUE} yields a posterior predictive interval
    for \eqn{\theta}. Using the \code{individual}
    argument, you can request individual effects'
    (\eqn{\theta_i}{theta[i]}) posterior distributions. May be an integer
    number (\code{1,...,k}) giving the index, or a character string giving
    the label.} 
  \item{cond.moment}{a \code{function(tau, predict=FALSE, individual=FALSE, simplify=TRUE)}
    returning conditional moments (mean and standard deviation) of
    \eqn{\mu} as a function of \eqn{\tau}. Using the argument
    \code{predict=TRUE} yields (conditional) posterior predictive
    moments for \eqn{\theta}. Using the \code{individual}
    argument, you can request individual effects'
    (\eqn{\theta_i}{theta[i]}) posterior distributions. May be an integer
    number (\code{1,...,k}) giving the index, or a character string giving
    the label.}
  \item{I2}{a \code{function(tau)} returning the \sQuote{relative} heterogeneity
    \eqn{I^2} as a function of \eqn{\tau}.}
  \item{summary}{a \code{matrix} listing some summary statistics, namely
    marginal posterior mode, median, mean, standard deviation
    and a (shortest) 95\% credibility intervals,
    of the marginal posterior distributions of \eqn{\tau} and \eqn{\mu},
    and of the posterior predictive distribution of \eqn{\theta}.}
  \item{interval.type}{the \code{interval.type} input argument
    specifying the type of interval to be returned by default.}
  \item{ML}{a \code{matrix} giving joint and marginal maximum-likelihood
    estimates of \eqn{(\tau,\mu)}.}
  \item{MAP}{a \code{matrix} giving joint and marginal
    maximum-a-posteriori estimates of \eqn{(\tau,\mu)}.}
  \item{theta}{a \code{matrix} giving the \sQuote{shrinkage estimates},
    i.e, summary statistics of the trial-specific means
    \eqn{\theta_i}{theta[i]}.} 
  \item{marginal.likelihood}{the marginal likelihood of the data (this 
    number is only computed if proper effect and heterogeneity priors 
    are specified).}
  \item{bayesfactor}{Bayes factors for the two hypotheses of
    \eqn{\tau=0} and \eqn{\mu=0}. These depend on the marginal
    likelihood and hence can only be computed if proper effect and
    heterogeneity priors are specified.}
  \item{support}{a \code{matrix} giving the \eqn{\tau} support points used
    internally in the grid approximation, along with their associated
    weights, conditional mean and standard deviation of \eqn{\mu},
    and the standard deviation of the (conditional) predictive
    distribution of \eqn{\theta}.}
  \item{delta, epsilon}{the \sQuote{\code{delta}} and \sQuote{\code{epsilon}}
    input parameter determining numerical accuracy.}
  \item{rel.tol.integrate, abs.tol.integrate}{the input parameters
    determining the numerical accuracy of the internally used
    \code{\link{integrate}()} function.}
  \item{call}{an object of class \code{call} giving the
    function call that generated the \code{bayesmeta} object.}
  \item{init.time}{the computation time (in seconds) used to generate
    the \code{bayesmeta} object.}
  }
\references{
  C. Roever et al. Bayesian random-effects meta-analysis. (in
  preparation), 2017.

  C. Roever, T. Friede.
  \href{http://doi.org/10.1080/10618600.2016.1276840}{Discrete
  approximation of a mixture distribution via restricted divergence}.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{26}(1):217-222, 2017.

  J.O. Berger, J. Deely.
  \href{http://doi.org/10.1080/01621459.1988.10478606}{A Bayesian
  approach to ranking and selection of related means with alternatives
  to analysis-of-variance methodology}. \emph{Journal of the American
  Statistical Association}, \bold{83}(402):364-373, 1988. 

  O. Bodnar, A. Link, C. Elster.
  \href{http://doi.org/10.1214/14-BA933}{Objective Bayesian inference
  for a generalized marginal random effects model}. \emph{Bayesian
  Analysis}, \bold{11}(1):25-45, 2016.

  A. Gelman, J.B. Carlin, H.S. Stern, D.B. Rubin. \emph{Bayesian data
  analysis}. Chapman & Hall / CRC, Boca Raton, 1997.
  
  A. Gelman. Prior distributions for variance parameters in hierarchical
  models. \emph{Bayesian Analysis}, \bold{1}(3):515-534, 2006.

  J. Hartung, G. Knapp, B.K. Sinha. \emph{Statistical meta-analysis with
  applications}. Wiley, Hoboken, 2008.

  L.V. Hedges, I. Olkin. \emph{Statistical methods for meta-analysis}.
  Academic Press, San Diego, 1985

  A.E. Kass, R.E. Raftery.
  \href{http://doi.org/10.2307/2291091}{Bayes factors}.
  \emph{Journal of the American Statistical Association},
  \bold{90}(430):773-795, 1995.

  D.V. Lindley.
  \href{http://doi.org/10.1093/biomet/44.1-2.187}{A statistical paradox}.
  \emph{Biometrika}, \bold{44}(1/2):187-192, 1957.
  
  B. Neuenschwander, G. Capkun-Niggli, M. Branson, D.J. Spiegelhalter.
  \href{http://doi.org/10.1177/1740774509356002}{Summarizing
  historical information on controls in clinical
  trials}. \emph{Trials}. \bold{7}(1):5-18, 2010.

  A. O'Hagan, L. Pericchi.
  \href{http://doi.org/10.1214/11-BJPS164}{Bayesian heavy-tailed
  models and conflict resolution: A review}. \emph{Brazilian Journal of
  Probability and Statistics}. \bold{26}(4):372-401, 2012.

  S. Shalloway. \href{http://doi.org/10.1214/14-BA883}{The
  evidentiary credible region}. \emph{Bayesian Analysis},
  \bold{9}(4):909-922, 2014. 

  D.J. Spiegelhalter, K.R. Abrams, J.P.Myles. \emph{Bayesian approaches
  to clinical trials and health-care evaluation}. Wiley & Sons, 2004.
  
  R. Tibshirani.
  \href{http://doi.org/10.1093/biomet/76.3.604}{Noninformative priors
  for one parameter of many}. \emph{Biometrika}, \bold{76}(3):604-608,
  1989.

  W. Viechtbauer. \href{http://www.jstatsoft.org/v36/i03/}{Conducting
  meta-analyses in R with the metafor package}. \emph{Journal of
  Statistical Software}, \bold{36}(3):1-48, 2010.
}
\author{
  Christian Roever \email{christian.roever@med.uni-goettingen.de}
}
\seealso{
  \code{\link{forestplot.bayesmeta}}, \code{\link{plot.bayesmeta}},
  \code{\link[metafor]{escalc}},
  \code{\link[compute.es:compute.es-package]{compute.es}}.
}
\examples{
########################################
# example data by Snedecor and Cochran:
data("SnedecorCochran")

\dontrun{
# analysis using improper uniform prior
# (may take a few seconds to compute!):
ma01 <- bayesmeta(y=SnedecorCochran[,"mean"], sigma=sqrt(SnedecorCochran[,"var"]),
                  label=SnedecorCochran[,"no"])

# analysis using an informative prior
# (normal for mu and half-Cauchy for tau (scale=10))
# (may take a few seconds to compute!):
ma02 <- bayesmeta(y=SnedecorCochran[,"mean"], sigma=sqrt(SnedecorCochran[,"var"]),
                  label=SnedecorCochran[,"no"],
                  mu.prior.mean=50, mu.prior.sd=50,
                  tau.prior=function(x){return(dhalfcauchy(x, scale=10))})

# show some summary statistics:
print(ma01)
summary(ma01)

# show some plots:
forestplot(ma01)
plot(ma01)

# compare resulting marginal densities;
# the effect parameter (mu):
mu <- seq(30, 80, le=200)
plot(mu, ma02$dposterior(mu=mu), type="l", lty="dashed",
     xlab=expression("effect "*mu),
     ylab=expression("marginal posterior density"),
     main="Snedecor/Cochran example")
lines(mu, ma01$dposterior(mu=mu), lty="solid")

# the heterogeneity parameter (tau):
tau <- seq(0, 50, le=200)
plot(tau, ma02$dposterior(tau=tau), type="l", lty="dashed",
     xlab=expression("heterogeneity "*tau),
     ylab=expression("marginal posterior density"),
     main="Snedecor/Cochran example")
lines(tau, ma01$dposterior(tau=tau), lty="solid")

# compute posterior median relative heterogeneity I-squared:
ma01$I2(tau=ma01$summary["median","tau"])

# determine 90 percent upper limits on the heterogeneity tau:
ma01$qposterior(tau=0.90)
ma02$qposterior(tau=0.90)
# determine shortest 90 percent credibility interval for tau:
ma01$post.interval(tau.level=0.9, method="shortest")}


#####################################
# example data by Sidik and Jonkman:
data("SidikJonkman2007")
# add log-odds-ratios and corresponding standard errors:
sj <- SidikJonkman2007
sj <- cbind(sj, "log.or"=log(sj[,"lihr.events"])-log(sj[,"lihr.cases"]-sj[,"lihr.events"])
                             -log(sj[,"oihr.events"])+log(sj[,"oihr.cases"]-sj[,"oihr.events"]),
                "log.or.se"=sqrt(1/sj[,"lihr.events"] + 1/(sj[,"lihr.cases"]-sj[,"lihr.events"])
                                 + 1/sj[,"oihr.events"] + 1/(sj[,"oihr.cases"]-sj[,"oihr.events"])))

\dontrun{
# analysis using weakly informative half-normal prior
# (may take a few seconds to compute!):
ma03a <- bayesmeta(y=sj[,"log.or"], sigma=sj[,"log.or.se"],
                   label=sj[,"id.sj"],
                   tau.prior=function(t){dhalfnormal(t,scale=1)})

# alternatively: may utilize "metafor" package's "escalc()" function
# to compute log-ORs and standard errors:
if (require("metafor")) {
  es <- escalc(measure="OR",
               ai=lihr.events, n1i=lihr.cases,
               ci=oihr.events, n2i=oihr.cases,
               slab=id, data=SidikJonkman2007)
} else print("Sorry, 'metafor' package not installed!")
# apply "bayesmeta()" function directly to "escalc" object:
ma03b <- bayesmeta(es, tau.prior=function(t){dhalfnormal(t,scale=1)})
# "ma03a" and "ma03b" should be identical:
print(ma03a$summary)
print(ma03b$summary)
# compare to metafor's (frequentist) random-effects meta-analysis:
rma03a <- rma.uni(es)
rma03b <- rma.uni(es, method="EB", knha=TRUE)
# compare mu estimates (estimate and confidence interval):
plot(ma03b, which=3)
abline(v=c(rma03a$b, rma03a$ci.lb, rma03a$ci.ub), col="red", lty=c(1,2,2))
abline(v=c(rma03b$b, rma03b$ci.lb, rma03b$ci.ub), col="green3", lty=c(1,2,2))
# compare tau estimates (estimate and confidence interval):
plot(ma03b, which=4)
abline(v=confint(rma03a)$random["tau",], col="red", lty=c(1,2,2))       
abline(v=confint(rma03b)$random["tau",], col="green3", lty=c(1,3,3))       

# show heterogeneity's posterior density:
plot(ma03a, which=4, main="Sidik/Jonkman example")

# show some numbers (mode, median and mean):
abline(v=ma03a$summary[c("mode","median","mean"),"tau"], col="blue")

# compare with Sidik and Jonkman's estimates:
sj.estimates <- sqrt(c("MM"  = 0.429,   # method of moments estimator
                       "VC"  = 0.841,   # variance component type estimator
                       "ML"  = 0.562,   # maximum likelihood estimator
                       "REML"= 0.598,   # restricted maximum likelihood estimator
                       "EB"  = 0.703,   # empirical Bayes estimator
                       "MV"  = 0.818,   # model error variance estimator
                       "MVvc"= 0.747))  # a variation of the MV estimator
abline(v=sj.estimates, col="red", lty="dashed")}


###########################
# example data by Cochran:
data("Cochran1954")

\dontrun{
# analysis using improper uniform prior
# (may take a few seconds to compute!):
ma04 <- bayesmeta(y=Cochran1954[,"mean"], sigma=sqrt(Cochran1954[,"se2"]),
                  label=Cochran1954[,"observer"])

# show joint posterior density:
plot(ma04, which=2, main="Cochran example")
# show (known) true parameter value:
abline(h=161)

# pick a point estimate for tau:
tau <- ma04$summary["median","tau"]
# highlight two point hypotheses (fixed vs. random effects):
abline(v=c(0, tau), col="orange", lty="dotted", lwd=2)

# show marginal posterior density:
plot(ma04, which=3)
abline(v=161)
# show the conditional distributions of the effect mu
# at two tau values corresponding to fixed and random effects models:
cm <- ma04$cond.moment(tau=c(0,tau))
mu <- seq(130,200, le=200)
lines(mu, dnorm(mu, mean=cm[1,"mean"], sd=cm[1,"sd"]), col="orange", lwd=2)
lines(mu, dnorm(mu, mean=cm[2,"mean"], sd=cm[2,"sd"]), col="orange", lwd=2)

# determine a range of tau values:
tau <- seq(0, ma04$qposterior(tau=0.99), length=100)
# compute conditional posterior moments:
cm.overall <- ma04$cond.moment(tau=tau)
# compute study-specific conditional posterior moments:
cm.indiv <- ma04$cond.moment(tau=tau, individual=TRUE)
# show forest plot along with conditional posterior means:
par(mfrow=c(1,2))
  plot(ma04, which=1, main="Cochran 1954 example")
  matplot(tau, cm.indiv[,"mean",], type="l", lty="solid", col=1:ma04$k,
          xlim=c(0,max(tau)*1.2), xlab=expression("heterogeneity "*tau),
          ylab=expression("(conditional) shrinkage estimate E["*
                           theta[i]*"|"*list(tau, y, sigma)*"]"))
  text(rep(max(tau)*1.01, ma04$k), cm.indiv[length(tau),"mean",],
       ma04$label, col=1:ma04$k, adj=c(0,0.5))
  lines(tau, cm.overall[,"mean"], lty="dashed", lwd=2)
  text(max(tau)*1.01, cm.overall[length(tau),"mean"],
       "overall", adj=c(0,0.5))
par(mfrow=c(1,1))

# show the individual effects' posterior distributions:
theta <- seq(120, 240, le=300)
plot(range(theta), c(0,0.1), type="n", xlab=expression(theta[i]), ylab="")
for (i in 1:ma04$k) {
  # draw estimate +/- uncertainty as a Gaussian:
  lines(theta, dnorm(theta, mean=ma04$y[i], sd=ma04$sigma[i]), col=i+1, lty="dotted")
  # draw effect's posterior distribution:
  lines(theta, ma04$dposterior(theta=theta, indiv=i), col=i+1, lty="solid")
}
abline(h=0)
legend(max(theta), 0.1, legend=ma04$label, col=(1:ma04$k)+1, pch=15, xjust=1, yjust=1)
}
}
\keyword{ models }
