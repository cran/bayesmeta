#
#    bayesmeta, an R package for Bayesian random-effects meta-analysis.
#    Copyright (C) 2021  Christian Roever
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


bayesmeta <- function(y,...)
{
  UseMethod("bayesmeta")
}


bayesmeta.default <- function(y, sigma, labels=names(y),
                              tau.prior="uniform",
                              mu.prior=c("mean"=NA,"sd"=NA),
                              mu.prior.mean=mu.prior[1],
                              mu.prior.sd  =mu.prior[2],
                              interval.type = c("shortest", "central"),
                              delta=0.01, epsilon=0.0001,
                              rel.tol.integrate=2^16*.Machine$double.eps,
                              abs.tol.integrate=0.0,
                              tol.uniroot=rel.tol.integrate,...)
{
  ptm <- proc.time()
  y      <- as.vector(y)
  sigma  <- as.vector(sigma)
  labels <- as.vector(labels)
  # some preliminary sanity checks:
  stopifnot(is.vector(y), is.vector(sigma),
            all(is.finite(y)), all(is.finite(sigma)),
            length(sigma)==length(y),
            all(sigma>=0), sum(sigma==0)<=1,
            length(mu.prior)==2,
            length(mu.prior.mean)==1, length(mu.prior.sd)==1,
            is.na(mu.prior.mean) || is.finite(mu.prior.mean),
            is.na(mu.prior.sd) || (is.finite(mu.prior.mean) && (mu.prior.sd>0)),
            ((is.na(mu.prior.mean) & is.na(mu.prior.sd))
             || (is.finite(mu.prior.mean) & is.finite(mu.prior.sd))),
            (is.function(tau.prior) | (is.character(tau.prior) && (length(tau.prior)==1))))
  interval.type <- match.arg(interval.type)
  stopifnot(length(interval.type)==1, is.element(interval.type, c("shortest","central")))

  zerosigma <- (sigma == 0.0)
  #if (any(zerosigma)) warning("one of the supplied 'sigma' elements is zero!")
  
  k <- length(y)
  if (is.null(labels))
    labels <- sprintf("%02d", 1:k)
  sigma2hat <- (k-1)*sum(1/sigma^2) / (sum(1/sigma^2)^2 - sum(1/sigma^4)) # Higgins/Thompson (2002), eqn. (9)

  maxratio <- max(sigma) / min(sigma)
  if ((!any(zerosigma)) && (maxratio > 1000))
    warning(paste0("Ratio of largest over smallest standard error (sigma) is ", sprintf("%.0f",maxratio), ". Extreme values may lead to computational problems."))

  tau.prior.proper <- NA
  if (is.character(tau.prior)) {
    tau.prior <- match.arg(tolower(tau.prior),
                           c("uniform","jeffreys","shrinkage","dumouchel","bergerdeely","conventional","i2","sqrt"))
    tau.prior <- c("uniform"="uniform", "jeffreys"="Jeffreys",
                   "shrinkage"="shrinkage", "dumouchel"="DuMouchel",
                   "bergerdeely"="BergerDeely", "conventional"="conventional", "i2"="I2", "sqrt"="sqrt")[tau.prior]
    stopifnot(is.element(tau.prior, c("uniform", "Jeffreys", "shrinkage", "DuMouchel", "BergerDeely", "conventional", "I2", "sqrt")))
    if (tau.prior=="uniform") {             # uniform prior on tau:
      pdens <- function(t){d<-rep(1,length(t)); d[t<0]<-0; return(d)}
      attr(pdens, "bayesmeta.label") <- "uniform(min=0, max=Inf)"
    } else if (tau.prior=="Jeffreys") {     # Jeffreys prior:
      pdens <- function(t){return(apply(matrix(t,ncol=1),1,function(x){return(sqrt(sum((x/(sigma^2+x^2))^2)))}))}
      attr(pdens, "bayesmeta.label") <- "Jeffreys prior"
    } else if (tau.prior=="BergerDeely") {  # Berger/Deely prior:
      pdens <- function(t){return(apply(matrix(t,ncol=1),1,
                                        function(x){return(exp(log(x)-sum(log(sigma^2+x^2))/k))}))}
      attr(pdens, "bayesmeta.label") <- "Berger/Deely prior"
    } else if (tau.prior=="conventional") { # conventional prior:
      pdens <- function(t){return(apply(matrix(t,ncol=1),1,
                                        function(x){return(exp(log(x)-(3/(2*k))*sum(log(sigma^2+x^2))))}))}
      attr(pdens, "bayesmeta.label") <- "conventional prior"
    } else if (tau.prior=="I2") {           # uniform on I^2:
      pdens <- function(t){return(apply(matrix(t,ncol=1),1,
                                        function(x){return(2 * exp(log(sigma2hat)+log(x) - 2*log(sigma2hat+x^2)))}))}
      attr(pdens, "bayesmeta.label") <- "uniform prior on I-squared"
    } else if (tau.prior=="sqrt") {         # sqrt-prior:
      pdens <- function(t){d <- rep(0,length(t)); d[t>=0] <- t[t>=0]^(-0.5); return(d)}
      attr(pdens, "bayesmeta.label") <- "uniform prior on sqrt(tau)"
    } else {
      # harmonic mean of squared standard errors:
      s02 <- k/sum(1/sigma^2)
      if (tau.prior=="shrinkage") {         # "uniform shrinkage" prior:
        pdens <- function(t){return(apply(matrix(t,ncol=1),1,function(x){return(2*x*s02/(s02+x^2)^2)}))}
        attr(pdens, "bayesmeta.label") <- "uniform shrinkage prior"
      } else if (tau.prior=="DuMouchel") {  # DuMouchel prior:
        s0 <- sqrt(s02)
        pdens <- function(t){return(apply(matrix(t,ncol=1),1,function(x){return(s0/(s0+x)^2)}))}
        attr(pdens, "bayesmeta.label") <- "DuMouchel prior"
      } else warning("could not make sense of 'tau.prior' argument")
    }
    if (is.element(tau.prior, c("uniform", "Jeffreys", "BergerDeely", "sqrt")))
      tau.prior.proper <- FALSE
    tau.prior <- pdens
    rm("pdens")
  }

  tau.prior.integral <- 1.0  # heterogeneity prior's normalizing constant

  dprior <- function(tau=NA, mu=NA, log=FALSE)
  # prior density (marginal or joint)
  {
    if (all(is.na(tau))) { # marginal density for mu:
      if (is.finite(mu.prior.mean))
        result <- dnorm(mu, mean=mu.prior.mean, sd=mu.prior.sd, log=log)
      else
        result <- rep(ifelse(log, 0, 1), length(mu))
    }
    else if (all(is.na(mu))) { # marginal density for tau:
      if (log) {
        if (is.element("log", names(formals(tau.prior))))
          result <- apply(matrix(tau,ncol=1), 1, tau.prior, log=TRUE) - log(tau.prior.integral)
        else
          result <- log(apply(matrix(tau,ncol=1), 1, tau.prior)) - log(tau.prior.integral)
      }
      else result <- apply(matrix(tau,ncol=1), 1, tau.prior) / tau.prior.integral
    }
    else if (is.finite(mu.prior.mean) & !all(is.na(mu))) {  # joint density:
      if (is.element("log", names(formals(tau.prior))))
        result <- (apply(matrix(tau,ncol=1), 1, tau.prior, log=TRUE) - log(tau.prior.integral)
                   + dnorm(mu, mean=mu.prior.mean, sd=mu.prior.sd, log=TRUE))
      else
        result <- (log(apply(matrix(tau,ncol=1), 1, tau.prior)) - log(tau.prior.integral)
                   + dnorm(mu, mean=mu.prior.mean, sd=mu.prior.sd, log=TRUE))
      if (!log) result <- exp(result)
    }
    else {  # joint density (uniform on mu):
      if (log) {
        if (is.element("log", names(formals(tau.prior))))
          result <- apply(matrix(tau,ncol=1), 1, tau.prior, log=TRUE) - log(tau.prior.integral)
        else
          result <- log(apply(matrix(tau,ncol=1), 1, tau.prior)) - log(tau.prior.integral)
      }
      else result <- apply(matrix(tau,ncol=1), 1, tau.prior) / tau.prior.integral
    }
    return(result)
  }

  # check whether tau prior is proper
  # (unless obvious from above specification)
  # by trying to integrate numerically:
  if (is.na(tau.prior.proper)) {
    prior.int <- integrate(function(t){return(dprior(tau=t))},
                           lower=0, upper=Inf,
                           rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate,
                           stop.on.error=FALSE)
    tau.prior.proper <- ((prior.int$message == "OK") && (prior.int$value > 0))
    # if integral is substantially different from 1.0, apply correction:
    if (tau.prior.proper && (abs(prior.int$value - 1.0) > prior.int$abs.error))
      tau.prior.integral <- prior.int$value 
    rm("prior.int")
  }
  
  likelihood <- function(tau=NA, mu=NA, log=FALSE)
  # likelihood function (marginal or joint)
  {
    if (all(is.na(mu)) & all(is.na(tau))) {
      warning("need to supply at least either 'mu' or 'tau'")
      return()
    }
    else if (all(is.na(tau))) { # return marginal likelihood (numerical):
      loglikeli <- function(taumu)
      {
        t2s2 <- taumu[1]^2 + sigma^2
        return(-(k/2)*log(2*pi) - 0.5*sum(log(t2s2)) - 0.5*sum((y-taumu[2])^2 / t2s2))
      }
      marglikeli <- function(mu)
      {
        integrand <- function(t){return(exp(loglikeli(c(t,mu)) + dprior(tau=t,log=TRUE)))}
        int <- integrate(Vectorize(integrand),
                         rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate,
                         lower=0, upper=Inf, stop.on.error=FALSE)
        return(ifelse(int$message == "OK", int$value, NA))
      }
      result <- apply(matrix(mu,ncol=1), 1, marglikeli)
      if (log) result <- log(result)
    } else if (all(is.na(mu))) {  # return marginal likelihood (analytical):
      logmarglikeli <- function(t)
      {
        if (is.na(mu.prior.mean)) { # uniform prior
          yTerm   <- y
          varTerm <- t^2 + sigma^2
        }
        else {                      # conjugate normal prior
          yTerm   <- c(mu.prior.mean, y)
          varTerm <- c(mu.prior.sd^2, t^2 + sigma^2)
        }
        if ((t==0) & any(zerosigma)) {
          conditionalmean <- y[zerosigma]
        } else {
          conditionalmean <- sum(yTerm/varTerm) / sum(1/varTerm)
        }
        return(-0.5*((length(yTerm)-1) * log(2*pi)
                     + sum(log(varTerm))
                     + sum((yTerm-conditionalmean)^2 / varTerm)
                     + log(sum(1/varTerm))))
      }
      result <- apply(matrix(tau, ncol=1), 1, logmarglikeli)
      result[tau<0] <- -Inf
      if (!log) result <- exp(result)
    }
    else {  # return joint likelihood:
      loglikeli <- function(taumu)
      {
        t2s2 <- taumu[1]^2 + sigma^2
        return(-(k/2)*log(2*pi) - 0.5*sum(log(t2s2)) - 0.5*sum((y-taumu[2])^2 / t2s2))
      }
      result <- apply(cbind(tau,mu), 1, loglikeli)
      result[tau<0] <- -Inf
      if (!log) result <- exp(result)
    }
    return(result)
  }
  
  dposterior <- function(tau=NA, mu=NA, theta=mu, log=FALSE, predict=FALSE, individual=FALSE)
  # posterior density (marginal or joint)
  {
    if (all(is.na(mu))) mu <- theta
    stopifnot(length(individual)==1)
    if (! (is.logical(individual) && (!individual))) {
      indiv.logi <- TRUE
      if (is.numeric(individual))   indiv.which <- which(is.element(1:k, individual))
      if (is.character(individual)) indiv.which <- which(is.element(labels, match.arg(individual,labels)))
      if (length(indiv.which)==0) warning("cannot make sense of 'individual' argument: empty subset.")
    }
    else indiv.logi <- FALSE
    if (predict & indiv.logi) warning("need to specify either 'predict' or 'individual' argument, but not both.")
    if (all(is.na(mu)) & all(is.na(tau))) {
      warning("need to supply at least either 'mu' or 'tau'")
      return()
    }
    else if (all(is.na(tau))) { # return marginal posterior (effect mu, numerical):
      if (all(is.na(support))) {
        warning("'support' not initialized.")
        result <- NA
      }
      else {
        result <- numeric(length(mu))
        if (predict)           # posterior predictive distribution
          for (i in 1:length(mu))
            result[i] <- sum(support[,"weight"]
                             * dnorm(mu[i], mean=support[,"mean"], sd=support[,"sd.pred"]))
        else if (indiv.logi) { # individual-effect distribution
          musigma <-   conditionalmoment(tau=support[,"tau"], individual=indiv.which)
          for (i in 1:length(mu))
            result[i] <- sum(support[,"weight"]
                             * dnorm(mu[i], mean=musigma[,"mean"], sd=musigma[,"sd"]))          
        }
        else                   # (marginal) posterior distribution
          for (i in 1:length(mu))
            result[i] <- sum(support[,"weight"]
                             * dnorm(mu[i], mean=support[,"mean"], sd=support[,"sd"]))
        if (log) result <- log(result)
      }
    }
    else {
      if (predict) warning("'predict' argument ignored!")
      if (indiv.logi) warning("'individual' argument ignored!")
      if (all(is.na(mu))) {  # return marginal posterior (heterogeneity tau, analytical):
        logmargpost <- function(t)
        {
          return(ifelse(t<0, -Inf, log(tau.prior(t)) + likelihood(mu=NA, tau=t, log=TRUE)))
        }
        result <- apply(matrix(tau, ncol=1), 1, logmargpost) - log(integral)
        result[is.nan(result)] <- -Inf
        if (!log) result <- exp(result)
      }
      else {  # return joint posterior:
        loglikeli <- function(taumu)
        {
          t2s2 <- taumu[1]^2 + sigma^2
          return(-(k/2)*log(2*pi) - 0.5*sum(log(t2s2)) - 0.5*sum((y-taumu[2])^2 / t2s2))
        }
        logpost <- function(taumu)
        {
          return(dprior(taumu[1], taumu[2], log=TRUE) + loglikeli(taumu) - log(integral))
        }
        result <- apply(cbind(tau,mu), 1, logpost)
        result[tau<0] <- -Inf
        if (!log) result <- exp(result)
      }
    }
    return(result)
  }

  # compute marginal posterior density's normalizing constant:
  integral <- 1.0
  integral <- integrate(function(x){dposterior(tau=x)}, lower=0, upper=Inf,
                        rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value
  if ((!is.finite(integral)) || (integral <= 0))
    warning("failed integrating marginal posterior (tau)")

  pposterior <- function(tau=NA, mu=NA, theta=mu, predict=FALSE, individual=FALSE)
  # posterior cumulative distribution function (CDF) of tau
  {
    if (all(is.na(mu))) mu <- theta
    stopifnot(length(individual)==1)
    if (! (is.logical(individual) && (!individual))) {
      indiv.logi <- TRUE
      if (is.numeric(individual))   indiv.which <- which(is.element(1:k, individual))
      if (is.character(individual)) indiv.which <- which(is.element(labels, match.arg(individual,labels)))
      if (length(indiv.which)==0) warning("cannot make sense of 'individual' argument: empty subset.")
    }
    else indiv.logi <- FALSE
    if (predict & indiv.logi) warning("need to specify either 'predict' or 'individual' argument, but not both.")
    if (all(is.na(mu))) { # numerical integration for tau
      if (predict) warning("'predict' argument ignored!")
      if (indiv.logi) warning("'individual' argument ignored!")
      cdf <- function(x)
      # marginal CDF of tau
      {
        if (x<=0) p <- 0
        else if (x==Inf) p <- 1
        else p <- integrate(dposterior, lower=0, upper=x,
                            rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value
        return(p)
      }
      result <- apply(matrix(tau, ncol=1), 1, cdf)
    }
    else if (all(is.na(tau))) { # grid approximation for mu
      if (all(is.na(support))) {
        warning("'support' not initialized.")
        result <- NA
      }
      else if (indiv.logi) { # individual-effect distribution
        musigma <-   conditionalmoment(tau=support[,"tau"], individual=indiv.which)
        cdf <- function(x)
        {
          p <- sum(pnorm(x, mean=musigma[,"mean"], sd=musigma[,"sd"]) * support[,"weight"])
          return(p)
        }
        result <- apply(matrix(mu, ncol=1), 1, cdf)
      }
      else { # posterior or posterior predictive CDF:
        cdf <- function(x)
        {
          p <- sum(pnorm(x, mean=support[,"mean"], sd=support[,ifelse(predict, "sd.pred", "sd")])
                   * support[,"weight"])
          return(p)
        }
        result <- apply(matrix(mu, ncol=1), 1, cdf)
      }
    }
    else {
      warning("need to supply EITHER tau OR mu, not both.")
      result <- NA
    }
    return(result)
  }

  qposterior <- function(tau.p=NA, mu.p=NA, theta.p=mu.p, predict=FALSE, individual=FALSE)
  # posterior quantile function of tau or mu
  {
    if (all(is.na(mu.p))) mu.p <- theta.p
    stopifnot(length(individual)==1)
    if (! (is.logical(individual) && (!individual))) {
      indiv.logi <- TRUE
      if (is.numeric(individual))   indiv.which <- which(is.element(1:k, individual))
      if (is.character(individual)) indiv.which <- which(is.element(labels, match.arg(individual,labels)))
      if (length(indiv.which)==0) warning("cannot make sense of 'individual' argument: empty subset.")
    }
    else indiv.logi <- FALSE
    if (all(is.na(mu.p))) { # compute tau quantile
      stopifnot(all(tau.p>=0), all(tau.p<=1))
      if (predict) warning("'predict' argument ignored!")
      if (indiv.logi) warning("'individual' argument ignored!")
      upper <- 1
      if (any((tau.p<1) & (tau.p>0))) {
        maxp <- max(tau.p[tau.p<1])
        while (pposterior(tau=upper) < maxp) upper <- upper * 2
      }
      qfun <- function(p)
      {
        stopifnot(p>=0, p<=1)
        if (p==0) quant <- 0
        else {
          if (p==1) quant <- Inf
          else
            quant <- uniroot(function(xx){return(pposterior(tau=xx)-p)},
                             interval=c(0,upper), tol=tol.uniroot)$root
        }
        return(quant)
      }
      result <- apply(matrix(tau.p,ncol=1),1,qfun)
    }
    else if (all(is.na(tau.p))) { # compute mu quantile
      if (indiv.logi && any(zerosigma) && (indiv.which == which(zerosigma))){
        result <- rep(NA, length(mu.p))
        result[is.finite(mu.p)] <- y[zerosigma]
      }
      stopifnot(all(mu.p>=0), all(mu.p<=1))
      if (any((mu.p<1) & (mu.p>0))) {
        minp <- min(c(mu.p[mu.p>0], 1-mu.p[mu.p<1]))
        # derive minimum/maximum based on variance and Chebychev inequality:
        if (predict) {
          lower <- sumstats["mean","theta"] - sqrt(1/minp)*sumstats["sd","theta"]*1.1
          upper <- sumstats["mean","theta"] + sqrt(1/minp)*sumstats["sd","theta"]*1.1
        }
        else if (indiv.logi) {
          lower <- y[indiv.which] - sqrt(1/minp)*sigma[indiv.which]*1.1
          upper <- y[indiv.which] + sqrt(1/minp)*sigma[indiv.which]*1.1
        }
        else {
          lower <- sumstats["mean","mu"] - sqrt(1/minp)*sumstats["sd","mu"]*1.1
          upper <- sumstats["mean","mu"] + sqrt(1/minp)*sumstats["sd","mu"]*1.1
        }
        while (pposterior(mu=lower,predict=predict,individual=individual) > minp)
            lower <- lower-sumstats["sd","mu"]
        while (pposterior(mu=upper,predict=predict,individual=individual) < (1-minp))
            upper <- upper+sumstats["sd","mu"]
      }
      qfun <- function(p)
      {
        stopifnot(p>=0, p<=1)
        if (p==0) quant <- -Inf
        else {
          if (p==1) quant <- Inf
          else
            quant <- uniroot(function(yy){return(pposterior(mu=yy,predict=predict,individual=individual)-p)},
                             interval=c(lower,upper), tol=tol.uniroot)$root
        }
        return(quant)
      }
      result <- apply(matrix(mu.p,ncol=1),1,qfun)
    }
    else {
      warning("need to supply EITHER tau.p OR mu.p, not both.")
      result <- NA
    }
    return(result)
  }

  rposterior <- function(n=1, predict=FALSE, individual=FALSE, tau.sample=TRUE)
  # posterior random number generation for tau and mu
  {
    stopifnot(n>0, n==round(n), length(individual)==1,
              !is.logical(individual) || !individual)
    if (tau.sample) {  # draw joint, bivariate (tau,mu) pairs:
      samp <- matrix(NA, nrow=n, ncol=2, dimnames=list(NULL,c("tau","mu")))
      if (is.numeric(individual) | is.character(individual))
          colnames(samp)[2] <- "theta"
      u <- runif(n=n)
      samp[,"tau"] <- apply(matrix(u,ncol=1), 1, function(x){return(qposterior(tau.p=x))})
      cond.sample <- function(t)
      {
        cm <- conditionalmoment(t, predict=predict, individual=individual)
        return(rnorm(n=1, mean=cm[1], sd=cm[2]))
      }
      samp[,2] <- apply(matrix(samp[,"tau"],ncol=1), 1, cond.sample)
    } else {           # draw marginal, univariate (mu or theta) numbers:
      samp <- rep(NA, n)
      if (!predict & (is.logical(individual) && (!individual)))
        meansd <- support[,c("mean","sd")]
      else
        meansd <- conditionalmoment(support[,"tau"], predict=predict, individual=individual)
      for (i in 1:n) {
        j <- sample(1:nrow(support), 1, prob=support[,"weight"])
        samp[i] <- rnorm(1, mean=meansd[j,"mean"], sd=meansd[j,"sd"])
      }          
    }
    return(samp)
  }

  post.interval <- function(tau.level=NA, mu.level=NA, theta.level=mu.level,
                            method=interval.type,
                            predict=FALSE, individual=FALSE)
  # determine credibility/prediction/shrinkage interval (tau, mu or theta)
  {
    if (is.na(mu.level)) mu.level <- theta.level
    stopifnot((is.finite(tau.level) & ((tau.level>0) & (tau.level<1)))
              | (is.finite(mu.level) & ((mu.level>0) & (mu.level<1))))
    stopifnot(length(individual)==1)
    method <- match.arg(method, c("shortest","central","evidentiary"))
    if (all(is.na(mu.level))) {      # CI for heterogeneity tau:
      if (predict) warning("'predict' argument ignored!")
      if (! (is.logical(individual) && (!individual))) warning("'individual' argument ignored!")
      if (method=="central")         # central interval
        result <- qposterior(tau.p=c((1-tau.level)/2, 1-(1-tau.level)/2))
      else if (method=="shortest") { # shortest interval
        intwidth <- function(left)
        {
          pleft <- pposterior(tau=left)
          right <- qposterior(tau=tau.level+pleft)
          return(right-left)
        }
        opti <- optimize(intwidth, lower=0, upper=qposterior(tau=1-tau.level))$minimum
        # catch marginal minimum:
        if (intwidth(0) < intwidth(opti))
          result <- c(0, qposterior(tau=tau.level))
        else
          result <- c(opti, qposterior(tau=tau.level+pposterior(tau=opti)))
      }
      else {                         # evidentiary interval (tau)
        expectedLogLikeli <- function(left)
        {
          pleft <- pposterior(tau=left)
          right <- qposterior(tau=tau.level+pleft)
          expect <- integrate(function(x){return(likelihood(tau=x,log=TRUE)*dposterior(tau=x))},
                              lower=left, upper=right,
                              rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)
          return(expect$value)
        }
        opti <- optimize(expectedLogLikeli, maximum=TRUE,
                         lower=0, upper=qposterior(tau=1-tau.level))$maximum
        # catch marginal minimum:
        if (expectedLogLikeli(0) > expectedLogLikeli(opti))
          result <- c(0, qposterior(tau=tau.level))
        else
          result <- c(opti, qposterior(tau=tau.level+pposterior(tau=opti)))          
      }
    }
    else if (all(is.na(tau.level))) {  # CI for effect mu:
      if (method=="central")           # central interval
        result <- qposterior(mu.p=c((1-mu.level)/2, 1-(1-mu.level)/2), predict=predict, individual=individual)
      else if (method=="shortest") {   # shortest interval
        intwidth <- function(left)
        {
          pleft <- pposterior(mu=left, predict=predict, individual=individual)
          right <- qposterior(mu=mu.level+pleft, predict=predict, individual=individual)
          return(right-left)
        }
        opti <- optimize(intwidth,
                         lower=qposterior(mu=(1-mu.level)/50, predict=predict, individual=individual),
                         upper=qposterior(mu=1-mu.level, predict=predict, individual=individual))$minimum
        result <- c(opti, qposterior(mu=mu.level+pposterior(mu=opti, predict=predict, individual=individual), predict=predict, individual=individual))
      }
      else {                           # evidentiary interval (mu)
        if (predict) {
          warning("evidentiary prediction intervals are not implemented!")
          result <- c(NA,NA)
        }
        else if ((length(individual)>1) || (!is.logical(individual)) || individual) {
          warning("evidentiary shrinkage intervals are not implemented!")
          result <- c(NA,NA)
        } else {
          expectedLogLikeli <- function(left)
          {
            pleft <- pposterior(mu=left)
            right <- qposterior(mu=mu.level+pleft)
            expect <- integrate(function(x){return(likelihood(mu=x)*dposterior(mu=x))},
                                lower=left, upper=right,
                                rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)
            return(expect$value)
          }
          opti <- optimize(expectedLogLikeli, maximum=TRUE,
                           lower=qposterior(mu=(1-mu.level)/100), upper=qposterior(mu=1-mu.level))$maximum
          result <- c(opti, qposterior(mu=mu.level+pposterior(mu=opti)))
        }
      }
    }
    else {
      warning("need to supply EITHER tau.level OR mu.level, not both.")
      result <- c(NA,NA)
    }
    attr(result, "interval.type") <- method
    return(result)
  }
  
  conditionalmoment <- function(tau, predict=FALSE, individual=FALSE, simplify=TRUE)
  # compute conditional posterior moments (mean, sd) of mu for given value(s) of tau
  # if (predict==TRUE), moments of the /posterior predictive distribution/ are returned.
  # if (individual==TRUE), moments of the conditional posterior /of the study-specific effects/
  # (theta[i]) are returned.
  {
    # interpret the "individual" argument, derive a "logical" and "numeric" component:
    if (all(is.logical(individual))) { #  (logical "individual" argument)
      indiv.logi  <- individual[1]
      if (indiv.logi)
        indiv.which <- 1:k
      else
        indiv.which <- NA          
    }
    else { #  (numerical/character "individual" argument)
      indiv.logi  <- TRUE
      if (all(is.numeric(individual))) # number(s) provided
        indiv.which <- which(is.element(1:k, individual))
      else if (all(is.character(individual))) # label(s) provided
        indiv.which <- which(is.element(labels, match.arg(individual, labels, several.ok=TRUE)))
      else warning("cannot make sense of 'individual' argument: funny format.")
      if (length(indiv.which)==0)
        warning("cannot make sense of 'individual' argument: empty subset.")
    }
    if (predict & indiv.logi) warning("need to specify either 'predict' or 'individual' argument, but not both.")
    cm <- function(t)
    {
      if ((t==0) && any(zerosigma)) {
        conditionalmean <- y[zerosigma]
        conditionalvar  <- 0.0
      } else {
        if (is.na(mu.prior.mean)) { # uniform prior
          yTerm   <- y
          varTerm <- t^2 + sigma^2
        }
        else {                      # conjugate normal prior
          yTerm   <- c(mu.prior.mean, y)
          varTerm <- c(mu.prior.sd^2, t^2 + sigma^2)
        }
        conditionalvar  <- 1 / sum(1/varTerm)
        conditionalmean <- sum(yTerm/varTerm) * conditionalvar
      }
      return(c("mean"=conditionalmean,"sd"=sqrt(conditionalvar)))
    }
    # compute conditional moments of mu
    # (or predictive, if requested):
    if (!indiv.logi) {
      result <- t(apply(matrix(tau,ncol=1),1,cm))
      if (predict) result[,"sd"] <- sqrt(result[,"sd"]^2 + tau^2)
    }
    # compute conditional moments of individual, study-specific means:
    else {
      result <- array(NA, dim=c(length(tau), 2, length(indiv.which)),
                      dimnames=list(NULL, c("mean","sd"), labels[indiv.which]))
      musigma <- t(apply(matrix(tau,ncol=1),1,cm))
      # loop over estimates:
      zero <- (tau==0)
      for (i in 1:length(indiv.which)) {
        if (any(!zero)) { # tau > 0
          if (sigma[indiv.which[i]] > 0.0) {
            result[!zero,"mean",i] <- (y[indiv.which[i]]/sigma[indiv.which[i]]^2 + musigma[!zero,"mean"]/tau[!zero]^2) / (1/sigma[indiv.which[i]]^2+1/tau[!zero]^2)
            result[!zero,"sd",i]   <- sqrt(1 / (1/sigma[indiv.which[i]]^2+1/tau[!zero]^2)
                                           + (musigma[!zero,"sd"] * (1/tau[!zero]^2) / (1/sigma[indiv.which[i]]^2+1/tau[!zero]^2))^2)
          } else {
            result[!zero,"mean",i] <- y[indiv.which[i]]
            result[!zero,"sd",i]   <- 0.0
          }
        }
        if (any(zero)) {  # tau = 0
          result[zero,"mean",i] <- musigma[zero,"mean"]
          result[zero,"sd",i]   <- musigma[zero,"sd"]
        }
      }
      # simplify array to matrix (if possible and desired):
      if ((dim(result)[3] == 1) && (simplify)) result <- result[,,1]
    }
    return(result)
  }

  ISquared <- function(tau)
  # compute heterogeneity measure "I-squared" as a function of tau
  {
    I2 <- rep(NA, length(tau))
    I2[tau>=0] <- tau[tau>=0]^2 / (tau[tau>=0]^2 + sigma2hat)
    return(I2)
  }

  discretize <- function(delta=0.01, epsilon=0.0001, alldivs=TRUE)
  # discretize parameter space in order to approximate marginal (mixture) distribution
  # via a finite number of mixture components; see also:
  #   C. Roever, T. Friede.
  #   Discrete approximation of a mixture distribution via restricted divergence.
  #   Journal of Computational and Graphical Statistics, 26(1): 217-222, 2017.
  #   http://doi.org/10.1080/10618600.2016.1276840
  {
    symKL <- function(mean1, sd1, mean2, sd2)
    # (general) symmetrized KL-divergence of two normal distributions
    {
      stopifnot(sd1>0, sd2>0)
      return((mean1-mean2)^2 * 0.5 * (1/sd1^2 + 1/sd2^2) + (sd1^2-sd2^2)^2 / (2*sd1^2*sd2^2))
    }
    divergence <- function(tau1, tau2)
    # (symmetrized) divergence b/w two conditionals specified through corresponding tau values
    {
      if (!alldivs) { # evaluate divergence based on (conditional) mu posterior only:
        cm <- conditionalmoment(tau=c(tau1, tau2))
        div <- symKL(cm[1,"mean"], cm[1,"sd"], cm[2,"mean"], cm[2,"sd"])
      }
      else { # evaluate divergence based also on (conditional) individual-study
             # and predictive mu posterior, and determine the maximum divergence:
        # determine array of ALL conditional moments:
        cm <- array(c(conditionalmoment(tau=c(tau1, tau2)),
                      conditionalmoment(tau=c(tau1, tau2), predict=TRUE),
                      #conditionalmoment(tau=c(tau1, tau2), individual=TRUE)),
                      conditionalmoment(tau=c(tau1, tau2), individual=(1:k)[!zerosigma])),
                    #dim=c(2, 2, k+2))
                    dim=c(2, 2, sum(!zerosigma)+2))
        # determine individual divergences and their maximum:
        #div <- rep(NA, k+2)
        div <- rep(NA, sum(!zerosigma)+2)
        #for (i in 1:(k+2))
        for (i in 1:(sum(!zerosigma)+2))
          div[i] <- symKL(cm[1,1,i], cm[1,2,i], cm[2,1,i], cm[2,2,i])
        div <- max(div)
      }
      return(div)
    }
    # determine range of interest / upper bound on tau:
    if (any(zerosigma)) maxtau <- qposterior(tau.p=1-min(c(epsilon/2, 1-epsilon/2)))
    else                maxtau <- qposterior(tau.p=1-min(c(epsilon, 1-epsilon)))
    # special procedure for 1st bin
    # (want zero to be first reference point UNLESS one of the sigma[i] is zero):
    if (any(zerosigma)) tau <- qposterior(tau.p=min(c(epsilon/2, 1-epsilon/2)))
    else                tau <- 0.0
    # search for upper bin margin:
    upper <- 1
    diverg <- divergence(tau, upper)
    while (diverg <= delta) {
      upper <- upper*2
      diverg <- divergence(tau, upper)
    }
    ur <- uniroot(function(t){return(divergence(tau, t)-delta)},
                  lower=tau, upper=upper,
                  f.lower=-delta, f.upper=diverg-delta,
                  tol=tol.uniroot)
    prob1 <- 0.0
    prob2 <- pposterior(tau=ur$root)
    # store result for 1st bin:
    result <- matrix(c(tau, prob2-prob1),
                     nrow=1, ncol=2,
                     dimnames=list(NULL,c("tau","weight")))
    tau <- ur$root
    # determine following bins (2,...):
    bin <- 2
    while ((tau <= maxtau) | (bin<=2)) {  # (at least 2 support points)
      result <- rbind(result, rep(0,2))
      # determine bin's reference point:
      diverg <- divergence(tau, upper)
      while (diverg <= delta) {
        upper <- upper*2
        diverg <- divergence(tau, upper)
      }
      ur <- uniroot(function(t){return(divergence(tau, t)-delta)},
                    lower=tau, upper=upper,
                    f.lower=-delta, f.upper=diverg-delta,
                    tol=tol.uniroot)
      tau <- ur$root
      result[bin,"tau"] <- tau
      # determine bin's upper bound:
      diverg <- divergence(tau, upper)
      while (diverg <= delta) {
        upper <- upper*2
        diverg <- divergence(tau, upper)
      }
      ur <- uniroot(function(t){return(divergence(tau, t)-delta)},
                    lower=tau, upper=upper,
                    f.lower=-delta, f.upper=diverg-delta,
                    tol=tol.uniroot)
      tau <- ur$root
      # determine bin's weight:
      prob1 <- prob2
      prob2 <- pposterior(tau=tau)
      result[bin,"weight"] <- max(c(0.0, prob2-prob1))
      # sanity check (to catch possible uniroot() failure):
      if (result[bin,"tau"] <= result[bin-1,"tau"]) {
        warning("DIRECT grid setup seems to have failed.")
        tau <- maxtau + 1.0
      }
      bin <- bin+1
    }
    # re-normalize weights (if necessary):
    sw <- sum(result[,"weight"])
    if (sw <= 0) {
      result[,"weight"] <- rep(1/nrow(result), nrow(result))
      sw <- 1.0
    }
    if (sw != 1.0)
      result[,"weight"] <- result[,"weight"] / sw
    return(result)
  }

  accelerate <- FALSE  # flag to speed up computations at the cost of a few estimates and a little accuracy

  # compute set of tau support points & corresponding weights:
  support <- discretize(delta=delta, epsilon=epsilon, alldivs=ifelse(accelerate, FALSE, TRUE))
  rm(list=c("discretize"))
  #if (nrow(support) <= 10)
  #  warning(paste0("Discretization of heterogeneity parameter space yielded very few (only ",nrow(support),") support points.\n",
  #                 "This *may* indicate numerical problems, e.g. due to a very narrow heterogeneity prior. Please double-check."))
  support <- cbind(support,
                   conditionalmoment(support[,"tau"]),
                   "sd.pred"=conditionalmoment(support[,"tau"],predict=TRUE)[,"sd"])
  
  # compute tau posterior's summary statistics:
  sumstats <- matrix(NA, nrow=6, ncol=3,
                     dimnames=list(c("mode", "median", "mean","sd", "95% lower", "95% upper"),
                                   c("tau","mu","theta")))
  expectation <- try(integrate(function(x)return(dposterior(x)*x), lower=0, upper=Inf,
                               rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value, silent=TRUE)
  if (class(expectation)=="try-error") {
    expectation <- NA
    variance <- NA
  }
  else {
    variance <- try(integrate(function(x)return(dposterior(x)*(x-expectation)^2), lower=0, upper=Inf,
                              rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value,
                    silent=TRUE)
    if (class(variance)=="try-error")
      variance <- NA
  }
  sumstats[c("mean","sd"),"tau"] <- c(expectation, sqrt(variance))
  sumstats["median","tau"] <- qposterior(0.5)
  sumstats[c("95% lower", "95% upper"),"tau"] <- post.interval(tau=0.95, method=ifelse(accelerate, "central", interval.type))
  if (!accelerate) {
    if (! is.finite(dposterior(tau=0))) {
      sumstats["mode","tau"] <- 0
    } else {
      maxi <- optimize(function(x){return(dposterior(tau=x))},
                       lower=0, upper=qposterior(tau=0.9), maximum=TRUE)
      if ((maxi$maximum <= 2 * .Machine$double.eps^0.25)
          && (dposterior(0) > maxi$objective)) maxi <- list("maximum"=0)
      sumstats["mode","tau"] <- maxi$maximum
      rm(list=c("maxi"))
    }
  }
  rm(list=c("expectation", "variance"))
  
  # compute mu posterior's summary statistics:
  if (all(is.finite(support))) {
    sumstats["mean","mu"] <- sum(support[,"mean"] * support[,"weight"])
    sumstats["sd","mu"] <- sqrt(sum(((support[,"mean"]-sumstats["mean","mu"])^2 + support[,"sd"]^2) * support[,"weight"]))
    sumstats["median","mu"] <- qposterior(mu=0.5)
    sumstats[c("95% lower", "95% upper"),"mu"] <- post.interval(mu=0.95, method=ifelse(accelerate, "central", interval.type))
    if (!accelerate) {
      maxi <- optimize(function(x){return(dposterior(mu=x))},
                       lower=qposterior(mu=0.1), upper=qposterior(mu=0.9), maximum=TRUE)
      sumstats["mode","mu"] <- maxi$maximum
      rm(list=c("maxi"))
    }
  }

  # compute mu posterior-predictive distribution's summary statistics:
  if (all(is.finite(support))) {
    sumstats["mean","theta"] <- sumstats["mean","mu"]
    sumstats["sd","theta"] <- sqrt(sum(((support[,"mean"]-sumstats["mean","theta"])^2 + support[,"sd.pred"]^2) * support[,"weight"]))
    if (!accelerate) {
      sumstats["median","theta"] <- qposterior(mu=0.5, predict=TRUE)
      sumstats[c("95% lower", "95% upper"),"theta"] <- post.interval(mu=0.95, predict=TRUE)
      maxi <- optimize(function(x){return(dposterior(mu=x, predict=TRUE))},
                       lower=qposterior(mu=0.1, predict=TRUE), upper=qposterior(mu=0.9, predict=TRUE), maximum=TRUE)
      sumstats["mode","theta"] <- maxi$maximum
      rm(list=c("maxi"))
    }
  }

  # compute joint & marginal maximum-likelihood (ML) estimates:
  ml.estimate <- rbind("joint"=c("tau"=NA, "mu"=NA), "marginal"=c("tau"=NA, "mu"=NA))
  map.estimate <- rbind("joint"=c("tau"=NA, "mu"=NA),
                        "marginal"=c("tau"=sumstats["mode","tau"],
                                     "mu"=sumstats["mode","mu"]))
  if (!accelerate) {
    # search for joint ML:
    opti <- optim(sumstats["median",c("tau","mu")],
                  function(x){return(-likelihood(tau=x[1], mu=x[2], log=TRUE))})
    # possibly catch optimum at (tau=0) parameter space margin:
    if (any(zerosigma)) {
      ml.estimate["joint",] <- opti$par
    } else {
      upper <- opti$par[2] + 1
      while (likelihood(tau=0, mu=upper+1, log=TRUE) > likelihood(tau=0, mu=upper, log=TRUE))
        upper <- upper + 1
      lower <- opti$par[2] - 1
      while (likelihood(tau=0, mu=lower-1, log=TRUE) > likelihood(tau=0, mu=lower, log=TRUE))
        lower <- lower - 1
      maxi <- optimize(function(x){return(likelihood(tau=0, mu=x, log=TRUE))},
                       lower=lower, upper=upper, maximum=TRUE)
      if (likelihood(tau=opti$par[1], mu=opti$par[2], log=TRUE) > maxi$objective)
        ml.estimate["joint",] <- opti$par
      else
        ml.estimate["joint",] <- c(0, maxi$maximum)
    }
    # marginal ML (mu):
    upper <- ml.estimate["joint","mu"] + 1
    while (likelihood(mu=upper+1) > likelihood(mu=upper))
      upper <- upper + 1
    lower <- ml.estimate["joint","mu"] - 1
    while (likelihood(mu=lower-1) > likelihood(mu=lower))
      lower <- lower - 1
    maxi <- optimize(function(x){return(likelihood(mu=x))},
                     lower=lower, upper=upper, maximum=TRUE)
    ml.estimate["marginal","mu"] <- maxi$maximum
    
    # marginal ML (tau):
    upper <- ifelse(ml.estimate["joint","tau"] > 0, 2*ml.estimate["joint","tau"], 1.0)
    while (likelihood(tau=2*upper) > likelihood(tau=upper))
      upper <- upper * 2
    maxi <- optimize(function(x){return(likelihood(tau=x))},
                     lower=0, upper=upper, maximum=TRUE)
    if (all(!zerosigma) && (likelihood(tau=0) > maxi$objective))
      ml.estimate["marginal","tau"] <- 0.0
    else
      ml.estimate["marginal","tau"] <- maxi$maximum
    
    # compute joint maximum-a-posteriori (MAP) estimate:
    opti <- optim(sumstats["median",c("tau","mu")],
                  function(x){return(-dposterior(tau=x[1], mu=x[2], log=TRUE))})
    map.value <- dposterior(tau=opti$par[1], mu=opti$par[2])
    map.estimate["joint",] <- opti$par
      
    # possibly catch optimum at (tau=0) parameter space margin:
    if (all(!zerosigma)) {
      if (is.finite(tau.prior(0))) {
        upper <- sumstats["95% upper", "mu"]
        while (dposterior(mu=upper+1, tau=0) > dposterior(mu=upper, tau=0))
          upper <- upper + 1
        lower <- sumstats["95% lower", "mu"]
        while (dposterior(mu=lower-1, tau=0) > dposterior(mu=lower, tau=0))
          lower <- lower - 1
        maxi <- optimize(function(x){return(dposterior(mu=x, tau=0))},
                         lower=lower, upper=upper, maximum=TRUE)
        if (maxi$objective > map.value) {
          map.estimate["joint",] <- c(0, maxi$maximum)
          map.value <- maxi$objective
        }
      } else {
        map.value <- Inf
        map.estimate["joint","tau"] <- 0.0
      }
    }
      
    # possibly catch diverging posterior density
    if (! is.finite(map.value)) {
      map.estimate["joint","mu"] <- conditionalmoment(tau=map.estimate["joint","tau"])[1,"mean"]
    }

    # clean up:
    rm(list=c("opti","maxi","lower","upper","map.value"))
  }
  # compute "shrinkage" estimates of theta[i]:
  shrink <- matrix(NA, nrow=8, ncol=k,
                   dimnames=list(c("y","sigma","mode", "median", "mean","sd", "95% lower", "95% upper"),
                                 labels))
  shrink["y",] <- y
  shrink["sigma",] <- sigma
  if (all(is.finite(support)) & (!accelerate)) {
    for (i in 1:k) {
      if (zerosigma[i]) {
        shrink[c("mean", "median", "mode", "95% lower", "95% upper"), i] <- y[i]
        shrink["sd", i] <- 0.0
      } else {
        musigma <- conditionalmoment(support[,"tau"], individual=i)
        shrink["mean",i] <- sum(musigma[,"mean"] * support[,"weight"])
        shrink["sd",i]   <- sqrt(sum(((musigma[,"mean"]-shrink["mean",i])^2 + musigma[,"sd"]^2) * support[,"weight"]))
        shrink["median",i] <- qposterior(mu=0.5, individual=i)
        shrink[c("95% lower", "95% upper"),i] <- post.interval(mu=0.95, individual=i)
        maxi <- optimize(function(x){return(dposterior(mu=x, individual=i))},
                         lower=qposterior(mu=0.1), upper=qposterior(mu=0.9), maximum=TRUE)
        shrink["mode",i] <- maxi$maximum
      }
    }
  }
  
  # compute marginal likelihood & Bayes factors:
  marglik <- NA_real_
  bayesfactor <- matrix(NA_real_, nrow=2, ncol=2, dimnames=list(c("actual", "minimum"), c("tau=0","mu=0")))
  if (tau.prior.proper) {
    bayesfactor["minimum","mu=0"]  <- likelihood(mu=0) / likelihood(mu=ml.estimate["marginal","mu"])
  }
  # check for proper effect prior:
  if (is.finite(mu.prior.mean) 
      && is.finite(mu.prior.sd)
      && (mu.prior.sd > 0)) {
    bayesfactor["minimum","tau=0"] <- likelihood(tau=0) / likelihood(tau=ml.estimate["marginal","tau"])
    # check for proper heterogeneity prior:
    if (tau.prior.proper) {
      marglik.int <- integrate(function(t){return(likelihood(tau=t)*dprior(tau=t))},
                               rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate,
                               lower=0, upper=Inf)
      if (marglik.int$message == "OK") {
        marglik <- marglik.int$value
        bayesfactor["actual", "tau=0"] <- likelihood(tau=0) / marglik
        bayesfactor["actual", "mu=0"]  <- likelihood(mu=0) / marglik
      } else {
        attr(marglik, "NA.reason") <- "failed computing marginal likelihood"
      }
      rm("marglik.int")
    } else {
      attr(marglik, "NA.reason") <- "improper heterogeneity prior"
    }
  } else {
    attr(marglik, "NA.reason") <- "improper effect prior"
  }

  # compute posterior mean inverse-variance weights:
  ivweights <- function(tau, idx)
  # inverse-variance weights as a function of tau
  {
    stopifnot(length(tau)==1, tau>=0)
    if (is.finite(mu.prior.mean) && is.finite(mu.prior.sd)) {
      weights <- 1 / c(tau^2 + sigma^2, mu.prior.sd^2)
      names(weights) <- c(labels, "prior mean")
    } else {
      weights <- 1 / (tau^2 + sigma^2)
      names(weights) <- labels
    }
    weights <- weights / sum(weights)
    return(weights[idx])
  }
  ivweights <- Vectorize(ivweights, vectorize.args="tau")
  
  # compute posterior means:
  if (is.finite(mu.prior.mean) && is.finite(mu.prior.sd)) {
    meanweights <- rep(NA_real_, k+1)
    names(meanweights) <- c(labels, "prior mean")
  } else {
    meanweights <- rep(NA_real_, k)
    names(meanweights) <- labels
  }
  for (i in 1:length(meanweights)) {
    meanweights[i] <- integrate(function(t){ivweights(tau=t, idx=i) * dposterior(tau=t)},
                                lower=0, upper=Inf,
                                rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value
  }

  # compute "shrinkage weights":
  shrinkweights <- function(tau, study.idx, shrink.idx)
  # "study.idx"  :  which study's weight
  # "shrink.idx" :  which shrinkage estimate
  {
    stopifnot(length(tau)==1, tau>=0)
    # compute inverse variance weights:
    if (is.finite(mu.prior.mean) && is.finite(mu.prior.sd)) {
      ivw <- 1 / c(tau^2 + sigma^2, mu.prior.sd^2)
      names(ivw) <- c(labels, "prior mean")
    } else {
      ivw <- 1 / (tau^2 + sigma^2)
      names(ivw) <- labels
    }
    ivw <- ivw / sum(ivw)
    if (any(sigma==0)) { # catch zero sigma case:
      s0 <- which(sigma==0)[1]
      for (i in 1:length(ivw)) ivw[i] <- 0.0
      ivw[s0] <- 1.0
    }
    # compute shrinkage weights:
    if (tau > 0)
      swt <- (sigma[shrink.idx]^(-2)) / (sigma[shrink.idx]^(-2) + tau^(-2))
    else
      swt <- 0.0
    indic <- rep(0.0, length(ivw))
    indic[shrink.idx] <- 1.0
    swt <- swt*indic + (1-swt)*ivw
    swt <- swt / sum(swt)
    if (any(sigma==0)) { # catch zero sigma case:
      s0 <- which(sigma==0)[1]
      swt <- rep(0.0, length(ivw))
      swt[s0] <- 1.0
    }
    return(swt[study.idx])
  }
  shrinkweights <- Vectorize(shrinkweights, vectorize.args="tau")
  shrinkageweights <- matrix(NA_real_, nrow=length(meanweights), ncol=k,
                             dimnames=list("study"=names(meanweights),
                                           "shrinkage estimate"=labels))
  for (i in 1:nrow(shrinkageweights)) {
    for (j in 1:ncol(shrinkageweights)) {
      shrinkageweights[i,j] <- integrate(function(t){shrinkweights(tau=t,i,j)*dposterior(tau=t)},
                                         lower=0.0, upper=Inf,
                                         rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value
    }
  }
  
  ptm <- proc.time()[1] - ptm[1]
  
  muprior <- c(mu.prior.mean, mu.prior.sd)
  names(muprior) <- c("mean","sd")
  
  # assemble eventual result to be returned:
  result <- list("y"                   = y,
                 "sigma"               = sigma,
                 "labels"              = labels,
                 "k"                   = k,
                 "tau.prior"           = tau.prior,
                 "mu.prior"            = muprior,
                 "dprior"              = dprior,
                 "tau.prior.proper"    = tau.prior.proper,
                 "likelihood"          = likelihood,
                 "dposterior"          = dposterior,
                 "pposterior"          = pposterior,
                 "qposterior"          = qposterior,
                 "rposterior"          = rposterior,
                 "post.interval"       = post.interval,
                 "cond.moment"         = conditionalmoment,
                 "I2"                  = ISquared,
                 "summary"             = sumstats,
                 "interval.type"       = interval.type,
                 "ML"                  = ml.estimate,
                 "MAP"                 = map.estimate,
                 "theta"               = shrink,
                 "weights"             = meanweights,
                 "weights.theta"       = shrinkageweights,
                 "marginal.likelihood" = marglik,
                 "bayesfactor"         = bayesfactor,
                 "support"             = support,
                 "delta"               = delta,
                 "epsilon"             = epsilon,
                 "rel.tol.integrate"   = rel.tol.integrate,
                 "abs.tol.integrate"   = abs.tol.integrate,
                 "tol.uniroot"         = tol.uniroot,
                 "call"                = match.call(expand.dots=FALSE),
                 "init.time"           = c("seconds"=unname(ptm)))
  class(result) <- "bayesmeta"
  return(result)
}


bayesmeta.escalc <- function(y, labels=NULL, ...)
# apply "bayesmeta()" to an "escalc" object
{
  attri <- attributes(y)
  if (all(is.element(c("yi.names", "vi.names"), names(attri)))) { # (for recent "metafor" versions)
    var.names <- c(attri$yi.names, attri$vi.names)
  } else if (is.element("var.names", names(attri))) {             # (for older "metafor" versions)
    var.names <- attri$var.names
  } else {
    stop(paste("Cannont extract \"yi.names\" and \"vi.names\" (or \"var.names\") attribute(s) from escalc object ",
               "(check use of the \"var.names\" option).", sep=""))
  }
  stopifnot(length(var.names)==2, all(is.character(var.names)))
  if (!all(is.element(var.names, names(y)))) {
    stop(paste("Cannont find columns \"",
               var.names[1],"\" and/or \"",
               var.names[2],"\" in escalc object ",
               "(check use of the \"var.names\" option).", sep=""))
  }
  if (is.null(labels)) {
    if (is.element("slab", names(attributes(y[,var.names[1]]))))
      labels <- as.character(attr(y[,var.names[1]], "slab"))
  }
  result <- bayesmeta.default(y=as.vector(y[,var.names[1]]),
                              sigma=sqrt(as.vector(y[,var.names[2]])),
                              labels=labels, ...)
  result$call <- match.call(expand.dots=FALSE)
  return(result)
}


print.bayesmeta <- function(x,...)
# print a short summary
{
  cat(" 'bayesmeta' object.\n")
  cat(paste("\n",x$k," estimates:\n", sep=""))
  if (length(x$labels)>10)
    cat(paste(c(x$labels[1:10], "..."), collapse=", "))
  else
    cat(paste(x$labels[1:length(x$labels)], collapse=", "))
  cat("\n")
  cat(paste0("\ntau prior (", ifelse(x$tau.prior.proper, "proper", "improper"), "):\n"))
  if (is.element("bayesmeta.label", names(attributes(x$tau.prior))))
    cat(paste(attributes(x$tau.prior)[["bayesmeta.label"]], "\n"))
  else
    print(x$tau.prior)
  mu.prior.proper <- (is.finite(x$mu.prior["mean"]) && is.finite(x$mu.prior["sd"]) && (x$mu.prior["sd"] > 0))
  cat(paste0("\nmu prior (", ifelse(mu.prior.proper, "proper", "improper"), "):\n"))
  if (mu.prior.proper)
    cat(paste("normal(mean=",x$mu.prior["mean"],", sd=",x$mu.prior["sd"],")\n",sep=""))
  else
    cat("uniform(min=-Inf, max=Inf)\n")
  cat("\nML and MAP estimates:\n")
  print(rbind("ML joint"=x$ML["joint",],
              "ML marginal"=x$ML["marginal",],
              "MAP joint"=x$MAP["joint",],
              "MAP marginal"=x$MAP["marginal",]))
  cat("\nmarginal posterior summary:\n")
  print(x$summary[,c("tau","mu")])
  if (x$interval.type=="shortest")
    cat("\n(quoted intervals are shortest credible intervals.)\n")
  else
    cat("\n(quoted intervals are central, equal-tailed credible intervals.)\n")      
  invisible(x)
}


summary.bayesmeta <- function(object,...)
# print a longer summary
{
  cat(" 'bayesmeta' object.\n")
  data <- matrix(c(object$y, object$sigma), ncol=2, dimnames=list(object$labels, c("y","sigma")))
  cat(paste("data (", object$k, " estimates):\n",sep=""))
  if (nrow(data)<=20)
    print(data)
  else {
    print(data[1:20,])
    cat(paste(" [...]  (truncated;", object$k, "estimates total)\n"))
  }
  cat(paste0("\ntau prior (", ifelse(object$tau.prior.proper, "proper", "improper"), "):\n"))
  if (is.element("bayesmeta.label", names(attributes(object$tau.prior))))
    cat(paste(attributes(object$tau.prior)[["bayesmeta.label"]], "\n"))
  else
    print(object$tau.prior)
  mu.prior.proper <- (is.finite(object$mu.prior["mean"]) && is.finite(object$mu.prior["sd"]) && (object$mu.prior["sd"] > 0))
  cat(paste0("\nmu prior (", ifelse(mu.prior.proper, "proper", "improper"), "):\n"))
  if (mu.prior.proper)
    cat(paste("normal(mean=",object$mu.prior["mean"],", sd=",object$mu.prior["sd"],")\n",sep=""))
  else
    cat("uniform(min=-Inf, max=Inf)\n")
  cat("\nML and MAP estimates:\n")
  print(rbind("ML joint"=object$ML["joint",],
              "ML marginal"=object$ML["marginal",],
              "MAP joint"=object$MAP["joint",],
              "MAP marginal"=object$MAP["marginal",]))
  cat("\nmarginal posterior summary:\n")
  print(object$summary)
  if (object$interval.type=="shortest")
    cat("\n(quoted intervals are shortest credible intervals.)\n")
  else
    cat("\n(quoted intervals are central, equal-tailed credible intervals.)\n")      
  if (any(is.finite(object$bayesfactor))) {
    cat("\nBayes factors:\n")
    print(object$bayesfactor)
  }
  cat("\nrelative heterogeneity I^2 (posterior median):", object$I2(tau=object$summary["median","tau"]), "\n")
  invisible(object)
}


plot.bayesmeta <- function(x, main=deparse(substitute(x)),
                           which=1:4, prior=FALSE,
                           forest.margin=8,
                           mulim=c(NA,NA), taulim=c(NA,NA),
                           violin=FALSE,...)
# generate forest plot and joint and marginal density plots.
{
  q975 <- qnorm(0.975)
  stopifnot(length(mulim)==2, length(taulim)<=2)

  if (all(is.finite(taulim))) {
    if ((length(taulim)==2) && (taulim[1]>=0) && (taulim[2]>taulim[1]))
      taurange <- taulim
    else if ((length(taulim)==1) && (taulim>0))
      taurange <- c(0, taulim)
    else
      taurange <- c(0, x$qposterior(tau=0.995)*1.1)
  } else {
    taurange <- c(0, x$qposterior(tau=0.995)*1.1)
  }

  if (all(is.finite(mulim)) && (mulim[1] < mulim[2])) {
    murange <- mulim
  } else {
    murange <- x$qposterior(mu=c(0.005,0.995))
    murange <- murange + c(-1,1)*diff(murange)*0.05
  }
  
  forestPlot <- function(x, main="", violin=violin, leftmargin=8)
  {
    original.margins <- par("mar")
    plotmargins <- original.margins
    plotmargins[2] <- leftmargin
    # determine x-axis range:
    xrange <- range(c(x$y-q975*x$sigma, x$y+q975*x$sigma,
                      x$summary[c("95% lower", "95% upper"),c("mu","theta")]))
    # empty plot:
    par("mar"=plotmargins)
    plot(xrange, c(-2, x$k)+c(-1,1)*0.5,
         type="n", axes=FALSE,
         xlab="", ylab="", main=main)
    mtext(side=1, line=par("mgp")[1], expression("effect"))
    # add horizontal line dividing data and estimates:
    abline(h=0, col="lightgrey")
    # add vertical "posterior median effect" line:
    abline(v=x$summary["median","mu"], col="black", lty="12")
    if (violin) {
      ticklength <- 0.45 # (here: maximum height of gaussian)
      maxdens <- c(dnorm(0, sd=x$sigma),
                   x$dposterior(mu=x$summary["mode","mu"]),
                   x$dposterior(theta=x$summary["mode","theta"], predict=TRUE))
      relmaxdens <- maxdens / max(maxdens)
      # density for central 95%
      narg1 <- seq(-q975, q975, le=51)
      ndens1 <- dnorm(narg1) / dnorm(0)
      # density for tails out to +/- 6 sigma
      narg2 <- seq(q975, 6.0, le=40)
      ndens2 <- dnorm(narg2) / dnorm(0)
      
      relsigma <- x$sigma / min(x$sigma)
      for (i in 1:x$k) { # loop over estimates
        # right tail:
        polygon(x$y[i] + c(narg2,rev(narg2))*x$sigma[i],
                x$k-(i-1) + c(ndens2,-rev(ndens2)) * relmaxdens[i] * ticklength,
                col="grey", border=NA)
        # left tail:
        polygon(x$y[i] - c(narg2,rev(narg2))*x$sigma[i],
                x$k-(i-1) + c(ndens2,-rev(ndens2)) * relmaxdens[i] * ticklength,
                col="grey", border=NA)
        # central chunk:
        polygon(x$y[i] + c(narg1,rev(narg1))*x$sigma[i],
                x$k-(i-1) + c(ndens1,-ndens1) * relmaxdens[i] * ticklength,
                col="black", border=NA)
        # central vertical line:
        lines(c(x$y[i], x$y[i]), x$k-(i-1)+ c(-1,1) * relmaxdens[i] * ticklength, col="grey")
      }
    }
    else {
      ticklength <- 0.3
      # draw horizontal lines for individual-study confidence intervals:
      matlines(rbind(x$y-q975*x$sigma, x$y+q975*x$sigma),
               rbind(x$k:1, x$k:1), col=grey(0.3), lty="solid")
      # draw vertical lines for individual-study effect estimates:
      matlines(rbind(x$y, x$y),
               rbind((x$k:1)+ticklength, (x$k:1)-ticklength), col="black", lty="solid")
    }
    if (violin) {
      # draw blob for mean effect estimate:
      ticklength <- 0.45 # (here: maximum height of gaussian)
      quant <- x$summary[c("95% lower","95% upper"),"mu"]
      quant <- c(quant[1]-2*(x$summary["median","mu"]-quant[1]),
                 quant,
                 quant[2]+2*(quant[2]-x$summary["median","mu"]))
      arg1 <- seq(quant[1], quant[2], le=50)  # left tail
      arg2 <- seq(quant[2], quant[3], le=51)  # central chunk
      arg3 <- seq(quant[3], quant[4], le=50)  # right tail
      dmode <- x$dposterior(mu=x$summary["mode","mu"])
      dens1 <- x$dposterior(mu=arg1) / dmode * relmaxdens[x$k+1]
      dens2 <- x$dposterior(mu=arg2) / dmode * relmaxdens[x$k+1]
      dens3 <- x$dposterior(mu=arg3) / dmode * relmaxdens[x$k+1]
      polygon(c(arg1, rev(arg1)), -1+c(dens1, -rev(dens1))*ticklength, col="grey", border=NA)
      polygon(c(arg3, rev(arg3)), -1+c(dens3, -rev(dens3))*ticklength, col="grey", border=NA)
      polygon(c(arg2, rev(arg2)), -1+c(dens2, -rev(dens2))*ticklength, col="black", border=NA)
      dm <- x$dposterior(mu=x$summary["median","mu"]) / dmode * relmaxdens[x$k+1]
      lines(rep(x$summary["median","mu"],2), -1+c(-1,1)*dm*ticklength, col="grey")
      
      # draw blob for prediction interval:
      quant <- x$summary[c("95% lower","95% upper"),"theta"]
      quant <- c(quant[1]-2*(x$summary["median","theta"]-quant[1]),
                 quant,
                 quant[2]+2*(quant[2]-x$summary["median","theta"]))
      arg1 <- seq(quant[1], quant[2], le=50)
      arg2 <- seq(quant[2], quant[3], le=51)
      arg3 <- seq(quant[3], quant[4], le=50)
      dmode <- x$dposterior(theta=x$summary["mode","theta"], predict=TRUE)
      dens1 <- x$dposterior(theta=arg1, predict=TRUE) / dmode * relmaxdens[x$k+2]
      dens2 <- x$dposterior(theta=arg2, predict=TRUE) / dmode * relmaxdens[x$k+2]
      dens3 <- x$dposterior(theta=arg3, predict=TRUE) / dmode * relmaxdens[x$k+2]
      polygon(c(arg1, rev(arg1)), -2+c(dens1, -rev(dens1))*ticklength, col="grey", border=NA)
      polygon(c(arg3, rev(arg3)), -2+c(dens3, -rev(dens3))*ticklength, col="grey", border=NA)
      polygon(c(arg2, rev(arg2)), -2+c(dens2, -rev(dens2))*ticklength, col="black", border=NA)
      dm <- x$dposterior(theta=x$summary["median","theta"], predict=TRUE) / dmode * relmaxdens[x$k+2]
      lines(rep(x$summary["median","theta"],2), -2+c(-1,1)*dm*ticklength, col="grey")
    }
    else {
      # draw diamond for mean effect estimate:
      ticklength <- 0.4
      polygon(x$summary[c("95% lower", "median", "95% upper", "median"),"mu"],
              rep(-1,4)+c(0,1,0,-1)*ticklength,
              border=NA, col=grey(0.4))
      # draw rectangle for effect prediction interval:
      ticklength <- 0.2
      polygon(x$summary[c("95% lower", "95% lower", "95% upper", "95% upper"),"theta"],
              rep(-2,4)+c(-1,1,1,-1)*ticklength,
              border=NA, col=grey(0.4))
    }
    # add axes & bounding box:
    axis(2, at=x$k:1, labels=x$labels, las=1)
    axis(2, at= c(-1,-2), labels=c(expression("effect "*mu), expression("prediction "*theta[italic(k)+1])), las=1)
    axis(1); box()
    par("mar"=original.margins)
    invisible()
  }
  
  jointdensity <- function(x, main="")
  {
    # range of tau values:
    tau <- seq(taurange[1], taurange[2],le=50)
    # range of mu values:
    mu <- seq(murange[1], murange[2], le=50)
    # grid of tau/mu value combinations:
    taumu <- expand.grid(tau,mu)
    # evaluate posterior density at grid points:
    post <- matrix(x$dposterior(tau=taumu[,1], mu=taumu[,2], log=TRUE),
                   nrow=length(tau), ncol=length(mu))
    # determine MAP value:
    map.value <- x$dposterior(x$MAP["joint",1], x$MAP["joint",2], log=TRUE)
    map.Inf <- !is.finite(map.value)
    if (map.Inf) {
      map.value <- max(post[is.finite(post)])
      warning("Non-finite posterior density, no contour lines drawn.", call.=FALSE)
    }
    # draw greyscale image:
    image(tau, mu, exp(post), axes=FALSE,
          col=grey((seq(1,0,le=128))^2),
          breaks=seq(0,exp(map.value),length=129),
          xlab="", ylab="", main=main, sub="(joint posterior density)", cex.sub=0.8)
    # add the blue lines for conditional mean & 95% confidence bounds:
    tau2 <- seq(taurange[1], taurange[2],le=200)
    cm <- x$cond.moment(tau=tau2)
    lines(tau2, cm[,"mean"], col="blue")  
    lines(tau2, cm[,"mean"]-q975*cm[,"sd"], col="blue", lty="dashed")  
    lines(tau2, cm[,"mean"]+q975*cm[,"sd"], col="blue", lty="dashed")
    # add green lines for marginal means & 95% confidence bounds:
    abline(v=x$summary[c("95% lower", "median", "95% upper"),"tau"],
           h=x$summary[c("95% lower", "median", "95% upper"),"mu"],
           col="green2", lty=c("16", "44", "16"))
    # draw ML estimate:
    points(x$ML["joint",1], x$ML["joint",2], col="magenta", pch=4)
    # draw MAP estimate:
    points(x$MAP["joint",1], x$MAP["joint",2], col="red", pch=3)
    if (!map.Inf) {
      # add contour lines:
      contour(tau, mu, post-map.value, add=TRUE, col="red",
              levels=-0.5*qchisq(p=c(0.5, 0.9, 0.95, 0.99), df=2),
              labels=paste(c(50, 90, 95, 99),"%",sep=""))
    }
    # add axes, bounding box, labels, ...
    axis(1); axis(2); box()
    mtext(side=1, line=par("mgp")[1], expression("heterogeneity "*tau))
    mtext(side=2, line=par("mgp")[1], expression("effect "*mu))
    invisible()
  }

  mumarginal <- function(x, main="", priorline=prior)
  {
    # range of mu values:
    mu <- seq(murange[1]-diff(murange)*0.05, murange[2]+diff(murange)*0.05, le=200)
    # corresponding posterior density:
    dens <- x$dposterior(mu=mu)
    # empty plot:
    plot(murange, c(0,max(dens,na.rm=TRUE)), type="n", axes=FALSE,
         xlab="", ylab="", main=main)
    # light grey shaded contour for density across whole range:
    polygon(c(min(mu), mu, max(mu)), c(0,dens,0), border=NA, col=grey(0.90))
    # dark grey shaded contour for density within 95% bounds:
    indi <- ((mu>=x$summary["95% lower","mu"]) & (mu<=x$summary["95% upper","mu"]))
    polygon(c(rep(x$summary["95% lower","mu"],2), mu[indi], rep(x$summary["95% upper","mu"],2)),
            c(0, x$dposterior(mu=x$summary["95% lower","mu"]),
              dens[indi], x$dposterior(mu=x$summary["95% upper","mu"]), 0),
            border=NA, col=grey(0.80))
    # vertical line for posterior median:
    lines(rep(x$summary["median","mu"],2),
          c(0,x$dposterior(mu=x$summary["median","mu"])), col=grey(0.6))
    # actual density line:
    lines(mu, dens, col="black")
    # x-axis:
    abline(h=0, col=grey(0.40))
    # prior density (if requested):
    if (priorline) lines(mu, x$dprior(mu=mu), col="black", lty="dashed")
    # add axes, labels, bounding box, ...
    mtext(side=1, line=par("mgp")[1], expression("effect "*mu))
    mtext(side=2, line=par("mgp")[2], expression("marginal posterior density"))
    axis(1); box()
    invisible()
  }

  taumarginal <- function(x, main="", priorline=prior)
  {
    # range of tau values:
    tau <- seq(max(c(0,taurange[1]-0.1*diff(taurange))),
                   taurange[2]+0.1*diff(taurange), le=200)
    # corresponding posterior density:
    dens <- x$dposterior(tau=tau)
    # empty plot:
    maxdens <- max(dens[is.finite(dens)],na.rm=TRUE)
    plot(c(taurange[1],taurange[2]), c(0,maxdens),         
         type="n", axes=FALSE, xlab="", ylab="", main=main)
    # "fix" diverging density:
    dens[!is.finite(dens)] <- 10*maxdens
    # light grey shaded contour for density across whole range:
    polygon(c(0,tau,max(tau)), c(0,dens,0), border=NA, col=grey(0.90))
    # dark grey shaded contour for density within 95% bounds:
    indi <- ((tau>=x$summary["95% lower","tau"]) & (tau<=x$summary["95% upper","tau"]))
    polygon(c(rep(x$summary["95% lower","tau"],2), tau[indi], rep(x$summary["95% upper","tau"],2)),
            c(0, min(c(x$dposterior(tau=x$summary["95% lower","tau"]), 10*maxdens)),
              dens[indi], x$dposterior(tau=x$summary["95% upper","tau"]), 0),
            border=NA, col=grey(0.80))
    # vertical line at posterior median:
    lines(rep(x$summary["median","tau"],2), c(0,x$dposterior(tau=x$summary["median","tau"])), col=grey(0.6))
    # actual density line:
    lines(tau, dens, col="black")
    # x-axis:
    abline(h=0, v=0, col=grey(0.40))
    # prior density (if requested):
    if (priorline) lines(tau, x$dprior(tau=tau), col="black", lty="dashed")
    # add axes, labels, bounding box, ...
    mtext(side=1, line=par("mgp")[1], expression("heterogeneity "*tau))
    mtext(side=2, line=par("mgp")[2], expression("marginal posterior density"))
    axis(1); box()
    invisible()
  }

  # main function:
  stopifnot(all(is.element(which, 1:4)),
            length(which)==length(unique(which)))  
  par.ask <- par("ask")
  for (i in 1:length(which)) {
    if (which[i]==1) forestPlot(x, main, violin=violin, leftmargin=forest.margin)
    else if (which[i]==2) jointdensity(x, main)
    else if (which[i]==3) mumarginal(x, main, priorline=prior)
    else if (which[i]==4) taumarginal(x, main, priorline=prior)
    if (i==1) {
      on.exit(par(ask=par.ask))
      par(ask=TRUE)
    }
  }
  par(ask=par.ask)
  invisible(x)
}


dlomax <- function(x, shape=1, scale=1, log=FALSE)
# probability density function for Lomax distribution
{
  stopifnot(length(shape)==1, length(scale)==1,
            all(shape>0), all(scale>0))
  result <- rep(-Inf, length(x))
  result[x>=0] <- log(shape)-log(scale)-(shape+1)*log(1+x[x>=0]/scale)
  if (!log) result <- exp(result)
  return(result)
}

plomax <- function(q, shape=1, scale=1)
# cumulative probability distribution function (CDF) of a Lomax distribution
{
  stopifnot(length(shape)==1, length(scale)==1,
            all(shape>0), all(scale>0))
  result <- rep(NA, length(q))
  result[q<=0] <- 0
  result[q>0] <- 1-(1+q[q>0]/scale)^(-shape)
  return(result)
}

qlomax <- function(p, shape=1, scale=1)
# quantile function (inverse CDF) of a Lomax distribution
{
  stopifnot(length(shape)==1, length(scale)==1,
            all(shape>0), all(scale>0))
  result <- rep(NA, length(p))
  result[p==1] <- Inf
  result[p==0] <- 0
  proper <- (is.finite(p) & (p>0) & (p<1))
  result[proper] <- ((1-p[proper])^(-1/shape) - 1) * scale
  return(result)
}

rlomax <- function(n, shape=1, scale=1)
# random number generation for a Lomax distribution
# (based on inversion method)
{
  stopifnot(length(shape)==1, length(scale)==1,
            all(shape>0), all(scale>0))
  u <- runif(n)
  result <- qlomax(u, scale=scale, shape=shape)
  return(result)
}

elomax <- function(shape=1, scale=1)
# expectation of a Lomax distribution
{
  stopifnot(length(shape)==1, length(scale)==1,
            all(shape>0), all(scale>0))
  if (shape > 1)
    expectation <- scale / (shape-1)
  else
    expectation <- Inf
  return(expectation)
}

vlomax <- function(shape=1, scale=1)
# variance of a Lomax distribution
{
  stopifnot(length(shape)==1, length(scale)==1,
            all(shape>0), all(scale>0))
  if (shape > 2)
    variance <- elomax(shape,scale)^2 * (shape/(shape-2))
  else
    variance <- Inf
  return(variance)
}


dhalfnormal <- function(x, scale=1, log=FALSE)
# probability density function of a half-normal distribution
{
  stopifnot(scale>0)
  result <- rep(-Inf, length(x))
  result[x>=0] <- log(2) + dnorm(x[x>=0], mean=0, sd=scale, log=TRUE)
  if (!log) result <- exp(result)
  return(result)  
}

phalfnormal <- function(q, scale=1)
# cumulative distribution function (CDF) of a half-normal distribution
{
  stopifnot(scale>0)
  result <- rep(0, length(q))
  result[q>0] <- 2.0 * (pnorm(q[q>0], mean=0, sd=scale) - 0.5)
  return(result)  
}

qhalfnormal <- function(p, scale=1)
# quantile function (inverse CDF) of a half-normal distribution
{
  stopifnot(scale>0)
  result <- rep(NA, length(p))
  proper <- (is.finite(p) & (p>=0) & (p<=1))
  result[proper] <- qnorm(p[proper]/2.0 + 0.5, mean=0, sd=scale)
  return(result)  
}

rhalfnormal <- function(n, scale=1)
# random number generation for half-normal distribution
{
  stopifnot(scale>0)
  return(abs(rnorm(n=n, mean=0, sd=scale)))  
}

ehalfnormal <- function(scale=1)
# expectation of a half-normal distribution
{
  stopifnot(scale>0)
  return(scale*sqrt(2/pi))  
}

vhalfnormal <- function(scale=1)
# variance of a half-normal distribution
{
  stopifnot(scale>0)
  return(scale^2*(1-2/pi))  
}


dhalft <- function(x, df, scale=1, log=FALSE)
# probability density function for half-Student-t distribution
{
  stopifnot(scale>0, df>0)
  result <- rep(-Inf, length(x))
  result[x>=0] <- log(2) - log(scale) + dt(x[x>=0]/scale, df=df, log=TRUE)
  if (!log) result <- exp(result)
  return(result)  
}

phalft <- function(q, df, scale=1)
# cumulative distribution function (CDF) of a half-Student-t distribution
{
  stopifnot(scale>0, df>0)
  result <- rep(0, length(q))
  result[q>0] <- 2.0 * (pt(q[q>0]/scale, df=df) - 0.5)
  return(result)  
}

qhalft <- function(p, df, scale=1)
# quantile function (inverse CDF) of a half-Student-t distribution
{
  stopifnot(scale>0, df>0)
  result <- rep(NA, length(p))
  proper <- (is.finite(p) & (p>=0) & (p<=1))
  result[proper] <- qt(p[proper]/2.0 + 0.5, df=df) * scale
  return(result)  
}

rhalft <- function(n, df, scale=1)
# random number generation for half-Student-t distribution
{
  stopifnot(scale>0, df>0)
  return(abs(rt(n=n, df=df))*scale)  
}

ehalft <- function(df, scale=1)
# expectation of a half-Student-t distribution
{
  stopifnot(scale>0, df>0)
  if (df>1)
    expectation <- exp(log(2)+log(scale)+0.5*log(df)-0.5*log(pi)+lgamma((df+1)/2)-lgamma(df/2)-log(df-1))
  else
    expectation <- Inf
  return(expectation)
}

vhalft <- function(df, scale=1)
# variance of a half-Student-t distribution
{
  stopifnot(scale>0, df>0)
  if (df>2)
    variance <- scale^2 * ((df / (df-2)) - ehalft(df=df, scale=1.0)^2)
  else
    variance <- Inf
  return(variance)
}


dhalfcauchy <- function(x, scale=1, log=FALSE)
# probability density function of a half-Cauchy distribution
{
  stopifnot(scale>0)
  result <- rep(-Inf, length(x))
  result[x>=0] <- log(2) + dcauchy(x[x>=0], location=0, scale=scale, log=TRUE)
  if (!log) result <- exp(result)
  return(result)  
}

phalfcauchy <- function(q, scale=1)
# cumulative distribution function (CDF) of a half-Cauchy distribution
{
  stopifnot(scale>0)
  result <- rep(0, length(q))
  result[q>0] <- 2.0 * (pcauchy(q[q>0], location=0, scale=scale) - 0.5)
  return(result)  
}

qhalfcauchy <- function(p, scale=1)
# quantile function (inverse CDF) of a half-Cauchy distribution
{
  stopifnot(scale>0)
  result <- rep(NA, length(p))
  proper <- (is.finite(p) & (p>=0) & (p<=1))
  result[proper] <- qcauchy(p[proper]/2.0 + 0.5, location=0, scale=scale)
  return(result)  
}

rhalfcauchy <- function(n, scale=1)
# random number generation for half-Cauchy distribution
{
  stopifnot(scale>0)
  return(abs(rcauchy(n=n, location=0, scale=scale)))
}

ehalfcauchy <- function(scale=1)
# expectation of a half-Cauchy distribution
{
  stopifnot(scale>0)
  return(Inf)
}

vhalfcauchy <- function(scale=1)
# variance of a half-Cauchy distribution
{
  stopifnot(scale>0)
  return(Inf)
}


dhalflogistic <- function(x, scale=1, log=FALSE)
# probability density function of a half-logistic distribution
{
  stopifnot(scale>0)
  result <- rep(-Inf, length(x))
  result[x>=0] <- log(2) + dlogis(x[x>=0], location=0, scale=scale, log=TRUE)
  if (!log) result <- exp(result)
  return(result)  
}

phalflogistic <- function(q, scale=1)
# cumulative distribution function (CDF) of a half-logistic distribution
{
  stopifnot(scale>0)
  result <- rep(0, length(q))
  result[q>0] <- 2.0 * (plogis(q[q>0], location=0, scale=scale) - 0.5)
  return(result)  
}

qhalflogistic <- function(p, scale=1)
# quantile function (inverse CDF) of a half-logistic distribution
{
  stopifnot(scale>0)
  result <- rep(NA, length(p))
  proper <- (is.finite(p) & (p>=0) & (p<=1))
  result[proper] <- qlogis(p[proper]/2.0 + 0.5, location=0, scale=scale)
  return(result)  
}

rhalflogistic <- function(n, scale=1)
# random number generation for half-logistic distribution
{
  stopifnot(scale>0)
  return(abs(rlogis(n=n, location=0, scale=scale)))
}

ehalflogistic <- function(scale=1)
# expectation of a half-logistic distribution
{
  stopifnot(scale>0)
  return(scale * log(4))
}

vhalflogistic <- function(scale=1)
# variance of a half-logistic distribution
{
  stopifnot(scale>0)
  return(scale^2 * (((pi^2)/3) - log(4)^2))
}


drayleigh <- function(x, scale=1, log=FALSE)
# probability density function of the Rayleigh distribution
{
  stopifnot(scale>0)
  result <- rep(-Inf, length(x))
  result[x>0] <- log(x[x>0])-2*log(scale)-0.5*(x[x>0]/scale)^2
  if (!log) result <- exp(result)
  return(result)
}

prayleigh <- function(q, scale=1)
# cumulative distribution function (CDF) of the Rayleigh distribution
{
  stopifnot(scale>0)
  result <- rep(0,length(q))
  result[q>0] <- 1-exp(-0.5*(q[q>0]/scale)^2)
  return(result)
}

qrayleigh <- function(p, scale=1)
# quantile function (inverse CDF) of the Rayleigh distribution
{
  stopifnot(scale>0)
  proper <- (is.finite(p) & (p>=0) & (p<=1))
  result <- rep(NA, length(p))
  result[proper] <- scale * sqrt(-2*log(1-p[proper]))
  return(result)
}

rrayleigh <- function(n, scale=1)
# random number generation for the Rayleigh distribution
{
  stopifnot(scale>0)
  return(sqrt(rexp(n, rate=1/(2*scale^2))))
}

erayleigh <- function(scale=1)
# expectation of a Rayleigh distribution
{
  stopifnot(scale>0)
  return(scale * sqrt(pi/2))
}

vrayleigh <- function(scale=1)
# variance of a Rayleigh distribution
{
  stopifnot(scale>0)
  return(scale^2 * (4-pi)/2)
}


dinvchi <- function(x, df, scale=1, log=FALSE)
# probability density function of a scaled inverse chi distribution
{
  stopifnot(length(df)==1, is.finite(df), df>0,
            length(scale)==1, is.finite(scale), scale>0)
  idx <- (is.finite(x) & (x>0))
  ldens <- rep(-Inf, length(x))
  ldens[idx] <- -(df/2-1)*log(2)-lgamma(df/2)-log(scale) -(df+1)*(log(x[idx])-log(scale)) - 0.5*(scale/x[idx])^2
  if (log) result <- ldens
  else     result <- exp(ldens)
  return(result)
}

pinvchi <- function(q, df, scale=1.0, lower.tail=TRUE, log.p=FALSE)
# cumulative distribution function (CDF) of a scaled inverse chi distribution
{
  stopifnot(length(df)==1, is.finite(df), df>0,
            length(scale)==1, is.finite(scale), scale>0)
  return(stats::pchisq(q=(q/scale)^(-2), df=df, lower.tail=(!lower.tail), log.p=log.p))
}

qinvchi <- function(p, df, scale=1.0, lower.tail=TRUE, log.p=FALSE)
# quantile function (inverse CDF) of a scaled inverse chi distribution
{
  stopifnot(length(df)==1, is.finite(df), df>0,
            length(scale)==1, is.finite(scale), scale>0)
  return(scale*stats::qchisq(p=p, df=df, lower.tail=(!lower.tail), log.p=log.p)^(-0.5))
}

rinvchi <- function(n=1, df, scale=1.0)
# random number generation for a scaled inverse chi distribution
{
  stopifnot(length(df)==1, is.finite(df), df>0,
            length(scale)==1, is.finite(scale), scale>0)
  return(scale * stats::rchisq(n=n, df=df)^(-0.5))
}

einvchi <- function(df, scale=1)
# expectation of a scaled inverse chi distribution
{
  stopifnot(length(scale)==1, is.finite(scale), scale>0,
            length(df)==1, is.finite(df), df>0)
  if (df<=1) result <- Inf
  else       result <- exp(log(scale)+lgamma((df-1)/2)-lgamma(df/2)-0.5*log(2))
  return(result)
}

vinvchi <- function(df, scale=1)
# variance of a scaled inverse chi distribution
{
  stopifnot(length(scale)==1, is.finite(scale), scale>0,
            length(df)==1, is.finite(df), df>0)
  if (df<=2) result <- Inf
  else       result <- scale^2/(df-2) - einvchi(df=df, scale=scale)^2
  return(result)
}


# Turner & al. prior data:
# ========================
# list of 5 possible intervention comparison types:
ctypes <- c("pharmacological vs. placebo / control",
            "pharmacological vs. pharmacological",
            "non-pharmacological vs. placebo / control",
            "non-pharmacological vs. pharmacological",
            "non-pharmacological vs. non-pharmacological")
# list of 16 possible outcome types:
otypes <- c("all-cause mortality",
            "obstetric outcomes",
            "cause-specific mortality / major morbidity event / composite (mortality or morbidity)",
            "resource use / hospital stay / process",
            "surgical / device related success / failure",
            "withdrawals / drop-outs",
            "internal / structure-related outcomes",
            "general physical health indicators",
            "adverse events",
            "infection / onset of new disease",
            "signs / symptoms reflecting continuation / end of condition",
            "pain",
            "quality of life / functioning (dichotomized)",
            "mental health indicators",
            "biological markers (dichotomized)",
            "subjective outcomes (various)")  
# matrix of mu values (see Tab.IV):
meanmat <- matrix(-c(3.95, 3.52, 3.71, 2.34, 2.14, 2.99, 2.71, 2.29, 1.87, 2.49, 2.06, 1.83, 2.54, 2.12, 1.77, 2.70,
                     4.18, 3.75, 3.95, 2.58, 2.37, 3.23, 2.94, 2.53, 2.10, 2.73, 2.29, 2.06, 2.78, 2.35, 2.00, 2.93,
                     4.17, 3.74, 3.93, 2.56, 2.36, 3.21, 2.93, 2.51, 2.10, 2.71, 2.28, 2.05, 2.77, 2.34, 1.99, 2.92,
                     2.92, 2.49, 2.68, 1.31, 1.11, 1.96, 1.67, 1.26, 0.84, 1.46, 1.03, 0.80, 1.51, 1.09, 0.74, 1.67,
                     3.50, 3.08, 3.27, 1.90, 1.69, 2.55, 2.26, 1.85, 1.43, 2.05, 1.61, 1.38, 2.10, 1.67, 1.33, 2.26),
                  nrow=16, ncol=5,
                 dimnames=list(otypes, ctypes))
# matrix of sigma values (see Tab.IV):
sdmat <- matrix(c(1.34, 1.74, 1.74, 1.74, 1.74, 1.74, 1.74, 1.53, 1.52, 1.52, 1.51, 1.52, 1.54, 1.53, 1.52, 1.52,
                  1.41, 1.79, 1.79, 1.79, 1.79, 1.79, 1.79, 1.58, 1.58, 1.58, 1.58, 1.58, 1.60, 1.60, 1.58, 1.58,
                  1.55, 1.91, 1.91, 1.91, 1.91, 1.91, 1.92, 1.72, 1.71, 1.71, 1.71, 1.71, 1.73, 1.72, 1.71, 1.71,
                  1.02, 1.50, 1.51, 1.50, 1.50, 1.51, 1.51, 1.25, 1.24, 1.24, 1.24, 1.25, 1.27, 1.27, 1.24, 1.25,
                  1.26, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.46, 1.45, 1.45, 1.45, 1.45, 1.47, 1.47, 1.45, 1.45),
                nrow=16, ncol=5,
                dimnames=list(otypes, ctypes))
# array of mu/sigma values:
TurnerEtAlParameters <- array(c(as.vector(meanmat), as.vector(sdmat)),
                              dim=c(16,5,2),
                              dimnames=list("outcome"=otypes, "comparison"=ctypes, "parameter"=c("mu","sigma")))
# remove obsolete objects:
rm(list=c("ctypes", "otypes", "meanmat", "sdmat"))


TurnerEtAlPrior <- function(outcome=c(NA,
                                      "all-cause mortality",
                                      "obstetric outcomes",
                                      "cause-specific mortality / major morbidity event / composite (mortality or morbidity)",
                                      "resource use / hospital stay / process",
                                      "surgical / device related success / failure",
                                      "withdrawals / drop-outs",
                                      "internal / structure-related outcomes",
                                      "general physical health indicators",
                                      "adverse events",
                                      "infection / onset of new disease",
                                      "signs / symptoms reflecting continuation / end of condition",
                                      "pain",
                                      "quality of life / functioning (dichotomized)",
                                      "mental health indicators",
                                      "biological markers (dichotomized)",
                                      "subjective outcomes (various)"),
                            comparator1=c("pharmacological", "non-pharmacological", "placebo / control"),
                            comparator2=c("pharmacological", "non-pharmacological", "placebo / control"))
#
# Function to return prior parameters, densities, etc. as proposed in
#
#   Turner et al.  Predictive distributions for between-study heterogeneity
#   and simple methods for their application in Bayesian meta-analysis.
#   Statisics in Medicine 4(6):984-998, 2015.
#
# (see Table IV).
#
{
  param <- matrix(NA, nrow=2, ncol=2,
                  dimnames=list(c("tau^2","tau"),
                                "log-normal"=c("mu","sigma")))
  if (all(!is.na(outcome))) {
    # match the provided arguments:
    outcome     <- match.arg(outcome)
    comparator1 <- match.arg(comparator1)
    comparator2 <- match.arg(comparator2)
    # list of 3 possible comparators:
    clist <- c("pharmacological", "non-pharmacological", "placebo / control")
    # index matrix to match pairs of comparators to one of five possible scenarios:
    cmatrix <- matrix(c(2,4,1,4,5,3,1,3,NA), nrow=3, ncol=3,
                      dimnames=list(clist, clist))
    # list of 5 possible intervention comparison types:
    ctypes <- dimnames(TurnerEtAlParameters)[["comparison"]]
    # figure out current comparison scenario:
    comparisontype <- ctypes[cmatrix[comparator1, comparator2]]
    # assemble function output:
    param["tau^2","mu"]    <- TurnerEtAlParameters[outcome, comparisontype, "mu"]
    param["tau^2","sigma"] <- TurnerEtAlParameters[outcome, comparisontype, "sigma"]
  } else {  # the "marginal" setting, see p.993
    outcome <- comparisontype <- "any"
    param <- matrix(NA, nrow=2, ncol=2,
                    dimnames=list(c("tau^2","tau"),
                                  "log-normal"=c("mu","sigma")))
    param["tau^2","mu"]    <- -2.56
    param["tau^2","sigma"] <- 1.74
  }
  param["tau","mu"]      <- param["tau^2","mu"] / 2
  param["tau","sigma"]   <- param["tau^2","sigma"] / 2
  result <- list("parameters"      = param,
                 "outcome.type"    = outcome,
                 "comparison.type" = comparisontype,
                 "dprior"          = function(tau, log=FALSE) {return(dlnorm(tau, meanlog=param["tau","mu"],
                                                                                  sdlog=param["tau","sigma"], log=log))},
                 "pprior"          = function(tau) {return(plnorm(tau, meanlog=param["tau","mu"], sdlog=param["tau","sigma"]))},
                 "qprior"          = function(p)   {return(qlnorm(p,   meanlog=param["tau","mu"], sdlog=param["tau","sigma"]))})
  attr(result$dprior, "bayesmeta.label") <- paste("log-normal(mu=",sprintf("%1.3f",param["tau","mu"]),
                                                  ", sigma=",sprintf("%1.3f",param["tau","sigma"]),")",sep="")
  return(result)
}



# Rhodes & al. prior data:
# ========================
# list of 3 possible intervention comparison types:
ctypes <- c("pharmacological vs. placebo / control",
            "pharmacological vs. pharmacological",
            "non-pharmacological (any)")

# list of 5 possible outcome types:
otypes <- c("obstetric outcome",
            "resource use and hospital stay / process",
            "internal and external structure-related outcome",            
            "general physical health and adverse event and pain and quality of life / functioning",
            "signs / symptoms reflecting continuation / end of condition and infection / onset of new acute / chronic disease",
            "mental health outcome",
            "biological marker",
            "various subjectively measured outcomes")

# matrix of location parameter values (see Tab.3):
locationmat <- matrix(-c(4.13, 2.55, 2.43, 3.16, 3.00, 2.99, 3.41, 2.76,
                         4.40, 2.83, 2.70, 3.44, 3.27, 3.27, 3.68, 3.03,
                         3.99, 2.41, 2.29, 3.02, 2.86, 3.85, 3.27, 2.62),
                      nrow=8, ncol=3,
                      dimnames=list(otypes, ctypes))

# matrix of location parameter values for RESPIRATORY DISEASES (see Tab.A.3.1):
locationmat.r <- matrix(-c(6.03, 4.46, 4.33, 5.07, 4.90, 4.90, 5.31, 4.66,
                           6.31, 4.73, 4.61, 5.34, 5.18, 5.17, 5.59, 4.94,
                           5.89, 4.32, 4.19, 4.93, 4.76, 4.76, 5.17, 4.52),
                        nrow=8, ncol=3,
                        dimnames=list(otypes, ctypes))

# matrix of location parameter values for CANCER (see Tab.A.3.2):
locationmat.c <- matrix(-c(1.57,-0.01, 0.13, 0.60, 0.44, 0.43, 0.85, 0.20,
                           1.85, 0.27, 0.14, 0.88, 0.71, 0.71, 1.13, 0.48,
                           1.43,-0.15,-0.27, 0.46, 0.30, 0.29, 0.71, 0.06),
                        nrow=8, ncol=3,
                        dimnames=list(otypes, ctypes))

# matrix of scale parameter values (see Tab.3):
scalemat <- matrix(c(2.34, 2.73, 2.50, 2.50, 2.50, 2.16, 2.83, 2.58,
                     2.31, 2.70, 2.46, 2.44, 2.47, 2.14, 2.78, 2.59,
                     2.11, 2.57, 2.32, 2.27, 2.33, 1.93, 2.66, 2.41),
                   nrow=8, ncol=3,
                   dimnames=list(otypes, ctypes))

# matrix of scale parameter values for RESPIRATORY DISEASES (see Tab.A.3.1):
scalemat.r <- matrix(c(2.36, 2.74, 2.51, 2.51, 2.50, 2.17, 2.83, 2.59,
                       2.31, 2.70, 2.46, 2.45, 2.47, 2.14, 2.78, 2.59,
                       2.21, 2.57, 2.33, 2.28, 2.33, 1.94, 2.66, 2.41),
                     nrow=8, ncol=3,
                     dimnames=list(otypes, ctypes))

# matrix of scale parameter values for CANCER (see Tab.A.3.2):
scalemat.c <- matrix(c(2.45, 2.83, 2.61, 2.61, 2.60, 2.28, 2.93, 2.68,
                       2.41, 2.79, 2.56, 2.55, 2.57, 2.25, 2.87, 2.68,
                       2.24, 2.68, 2.45, 2.40, 2.46, 2.08, 2.78, 2.53),
                     nrow=8, ncol=3,
                     dimnames=list(otypes, ctypes))

# array of mu/sigma values:
RhodesEtAlParameters <- array(c(as.vector(locationmat), as.vector(scalemat),
                                as.vector(locationmat.r), as.vector(scalemat.r),
                                as.vector(locationmat.c), as.vector(scalemat.c)),
                              dim=c(8,3,2,3),
                              dimnames=list("outcome"=otypes, "comparison"=ctypes,
                                            "parameter"=c("location","scale"),
                                            "medical area"=c("other","respiratory","cancer")))
# remove obsolete objects:
rm(list=c("ctypes", "otypes", "locationmat", "scalemat", "locationmat.r", "scalemat.r", "locationmat.c", "scalemat.c"))

RhodesEtAlPrior <- function(outcome=c(NA,
                                      "obstetric outcome",
                                      "resource use and hospital stay / process",
                                      "internal and external structure-related outcome",            
                                      "general physical health and adverse event and pain and quality of life / functioning",
                                      paste("signs / symptoms reflecting continuation / end of condition and infection",
                                            "/ onset of new acute / chronic disease"),
                                      "mental health outcome",
                                      "biological marker",
                                      "various subjectively measured outcomes"),
                            comparator1=c("pharmacological", "non-pharmacological", "placebo / control"),
                            comparator2=c("pharmacological", "non-pharmacological", "placebo / control"),
                            area=c("other","respiratory","cancer"))
#
# Function to return prior parameters as proposed in
#
#   Rhodes et al.
#   Predictive distributions were developed for the extent of heterogeneity in meta-analyses of continuous outcome data.
#   Journal of Clinical Epidemiology 68(1):52-60, 2015.
#
# (see Table 3).
#
# Technically, this function mostly retrieves the parameters from a pre-defined array,
# the 3-dimensional "RhodesEtAlParameters" array. (This is more convenient than having
# to deal with the array itself, mostly due to partial argument matching, etc.)
#
{
  # match the provided arguments:
  if (all(!is.na(outcome))) { # settings from Tables 3, A.3.1 or A.3.2
    outcome     <- match.arg(outcome)
    comparator1 <- match.arg(comparator1)
    comparator2 <- match.arg(comparator2)
    area        <- match.arg(area)
    # list of 3 possible comparators:
    clist <- c("pharmacological", "non-pharmacological", "placebo / control")
    # index matrix to match pairs of comparators to one of 3 possible scenarios:
    cmatrix <- matrix(c(2,3,1,3,3,3,1,3,NA), nrow=3, ncol=3,
                      dimnames=list(clist, clist))
    # list of 3 possible intervention comparison types:
    ctypes <- dimnames(RhodesEtAlParameters)[["comparison"]]
    # figure out current comparison scenario:
    comparisontype <- ctypes[cmatrix[comparator1, comparator2]]
    # assemble function output:
    location <- RhodesEtAlParameters[outcome, comparisontype, "location", area]
    scale    <- RhodesEtAlParameters[outcome, comparisontype, "scale", area]
  } else { # the general (marginal) setting; see beginning of Sec. 3.3
    outcome <- comparisontype <- area <- "any"
    location <- -3.44
    scale <- 2.59
  }
  param <- matrix(NA, nrow=2, ncol=2,
                  dimnames=list(c("tau^2","tau"),
                                "log-t(5)"=c("location","scale")))
  param["tau^2","location"] <- location
  param["tau^2","scale"]    <- scale
  param["tau","location"]   <- location / 2
  param["tau","scale"]      <- scale / 2
  rm(list=c("location","scale"))
  dlt5 <- function(x, location=0, scale=1, log=FALSE)
  # density of log-t distribution with 5 d.f.
  {
    stopifnot(is.finite(location), is.finite(scale), scale>0)
    logdensity <- rep(-Inf, length(x))
    proper <- (is.finite(x) & (x>0))
    logdensity[proper] <- dt((log(x[proper])-location)/scale, df=5, log=TRUE) - log(scale) - log(x[proper])
    if (log) return(logdensity)
    else     return(exp(logdensity))
  }
  plt5 <- function(x, location=0, scale=1)
  # cumulative distribution function (CDF) of log-t distribution with 5 d.f.
  {
    stopifnot(is.finite(location), is.finite(scale), scale>0)
    return(pt((log(x)-location)/scale, df=5))
  }
  qlt5 <- function(p, location=0, scale=1)
  # quantile function (inverse CDF) of log-t distribution with 5 d.f.
  {
    stopifnot(is.finite(location), is.finite(scale), scale>0)
    return(exp((qt(p, df=5)*scale)+location))
  }
  result <- list("parameters"      = param,
                 "outcome.type"    = outcome,
                 "comparison.type" = comparisontype,
                 "medical.area"    = area,
                 "dprior"          = function(tau, log=FALSE) {return(dlt5(tau, location=param["tau","location"],
                                                                                scale=param["tau","scale"], log=log))},
                 "pprior"          = function(tau) {return(plt5(tau, location=param["tau","location"], scale=param["tau","scale"]))},
                 "qprior"          = function(p)   {return(qlt5(p, location=param["tau","location"], scale=param["tau","scale"]))})
  attr(result$dprior, "bayesmeta.label") <- paste("log-Student-t(location=",sprintf("%1.3f",param["tau","location"]),
                                                  ", scale=",sprintf("%1.3f",param["tau","scale"]),", d.f.=5)",sep="")
  return(result)
}


forest.bayesmeta <- function(x, xlab="effect size", refline=0, cex=1,...)
# forest plot for a "bayesmeta" object
# based on the "metafor" package's plotting functions
{
  if (!requireNamespace("metafor", quietly=TRUE))
    stop("required 'metafor' package not available!")
  metafor::forest.default(x=x$y, sei=x$sigma,
                          showweight=FALSE,  # (IV-weights don't make sense here)
                          ylim=c(-x$k-4, 1),
                          level=95,          # (95% level is intentionally hard-coded)
                          refline=refline,
                          xlab=xlab,
                          slab=x$labels,
                          rows=seq(-2, -x$k - 1, by = -1),
                          cex=cex, ...)
  metafor::addpoly(x$summary["median","mu"], ci.lb=x$summary["95% lower","mu"], ci.ub=x$summary["95% upper","mu"],
                   rows = -x$k-2.5, mlab=expression("mean effect ("*mu*")"), level=95, cex=cex, ...)
  metafor::addpoly(x$summary["median","theta"], ci.lb=x$summary["95% lower","theta"], ci.ub=x$summary["95% upper","theta"],
                   rows = -x$k-3.5, mlab=expression("prediction ("*vartheta[k+1]*")"), level=95, cex=cex, ...)
  
  plotdata <- cbind("95% lower"=x$y-qnorm(0.975)*x$sigma, "estimate"=x$y, "95% upper"=x$y+qnorm(0.975)*x$sigma)
  plotdata <- rbind(plotdata, t(x$summary[c("95% lower","median","95% upper"),c("mu","theta")]))
  rownames(plotdata) <- c(x$labels, c("mean effect(mu)", "prediction (theta)"))
  invisible(plotdata)
}


forestplot.bayesmeta <- function(x, labeltext,
                                 exponentiate  = FALSE,
                                 prediction    = TRUE,
                                 shrinkage     = TRUE,
                                 heterogeneity = TRUE,
                                 digits        = 2,
                                 plot          = TRUE,
                                 fn.ci_norm, fn.ci_sum, col, legend, boxsize, ...)
#
# ARGUMENTS:
#   x            :  a "bayesmeta" object.
#   labeltext    :  you may provide an alternative "labeltext" argument here
#                   (see also the "forestplot()" help).
#   exponentiate :  flag indicating whether to exponentiate numbers (figure and table).
#   prediction   :  flag indicating whether to show prediction interval.
#   shrinkage    :  flag indicating whether to show shrinkage estimates.
#   digits       :  number of significant digits to be shown (based on standard deviations).
#   plot         :  flag you can use to suppress actual plotting.
#   ...          :  further arguments passed to the "forestplot" function.
#   
# VALUE:
#   a list with components
#     $data       :  the meta-analyzed data, and mean and prediction estimates
#     $shrinkage  :  the shrinkage estimates for each study
#     $labeltext  :  the "forestplot()" function's "labeltext" argument used internally
#     $forestplot :  the "forestplot()" function's returned value
#
{
  if (!requireNamespace("forestplot", quietly=TRUE))
    stop("required 'forestplot' package not available!")
  if (utils::packageVersion("forestplot") < "1.5.2")
    warning("you may need to update 'forestplot' to a more recent version (>=1.5.2).")
  # auxiliary function:
  decplaces <- function(x, signifdigits=3)
  # number of decimal places (after decimal point)
  # to be displayed if you want at least "signifdigits" digits
  {
    return(max(c(0, -(floor(log10(x))-(signifdigits-1)))))
  }
  # some sanity checks for the provided arguments:
  stopifnot(is.element("bayesmeta", class(x)),
            length(digits)==1, digits==round(digits), digits>=0,
            length(exponentiate)==1, is.logical(exponentiate),
            length(prediction)==1, is.logical(prediction),
            length(shrinkage)==1, is.logical(shrinkage),
            length(plot)==1, is.logical(plot))
  # plotting data (1) -- the quoted estimates:
  q95 <- qnorm(0.975)
  ma.dat <- rbind(NA,
                  cbind(x$y, x$y - q95*x$sigma, x$y + q95*x$sigma),
                  x$summary[c("median", "95% lower", "95% upper"),"mu"],
                  x$summary[c("median", "95% lower", "95% upper"),"theta"])
  colnames(ma.dat) <- c("estimate", "lower", "upper")
  rownames(ma.dat) <- c("", x$label, "mean", "prediction")
  if (! prediction) ma.dat <- ma.dat[-(x$k+3),]
  # plotting data (2) -- the shrinkage estimates:
  ma.shrink <- rbind(NA,
                     t(x$theta)[,c("median","95% lower","95% upper")],
                     NA,NA)
  colnames(ma.shrink) <- c("estimate", "lower", "upper")
  rownames(ma.shrink) <- c("", x$label, "mean", "prediction")
  if (! prediction) ma.shrink <- ma.shrink[-(x$k+3),]
  if (exponentiate) {
    ma.dat    <- exp(ma.dat)
    ma.shrink <- exp(ma.shrink)
  }
  # generate "labeltext" data table for plot (unless already provided):
  if (missing(labeltext)) {
    # determine numbers of digits based on standard deviations:
    if (exponentiate) {
      stdevs <- c(exp(x$y)*x$sigma,
                  exp(x$summary["median","mu"])*x$summary["sd","mu"])
      if (prediction) stdevs <- c(stdevs, exp(x$summary["median","theta"])*x$summary["sd","theta"])
    } else {
      stdevs <- c(x$sigma, x$summary["sd","mu"])
      if (prediction) stdevs <- c(stdevs, x$summary["sd","theta"])
    }
    stdevs <- abs(stdevs[is.finite(stdevs) & (stdevs != 0)])
    formatstring <- paste0("%.", decplaces(stdevs, digits), "f")
    # fill data table:
    labeltext <- matrix(NA_character_, nrow=nrow(ma.dat), ncol=3)
    labeltext[1,] <- c("study","estimate", "95% CI")
    labeltext[,1] <- c("study", x$labels, "mean", "prediction")[1:nrow(ma.dat)]
    for (i in 2:(nrow(ma.dat))) {
      labeltext[i,2] <- sprintf(formatstring, ma.dat[i,"estimate"])
      labeltext[i,3] <- paste0("[", sprintf(formatstring, ma.dat[i,"lower"]),
                                 ", ", sprintf(formatstring, ma.dat[i,"upper"]), "]")
    }
  }
  # add horizontal lines to plot:
  horizl <- list(grid::gpar(col="grey"), grid::gpar(col="grey"))
  names(horizl) <- as.character(c(2,x$k+2))
  # specify function(s) for drawing estimates / shrinkage estimates:
  if (missing(fn.ci_norm)) {
    if (shrinkage) {
      fn.ci_norm <- list(function(...) {forestplot::fpDrawPointCI(pch=15,...)},
                         function(...) {forestplot::fpDrawPointCI(pch=18,...)})
    } else {
      fn.ci_norm <- function(...) {forestplot::fpDrawPointCI(pch=15,...)}
    }
  }
  # specify function(s) for drawing summaries (diamond / bar):
  if (missing(fn.ci_sum)) {
    fn.ci_sum <- list(NULL)
    for (i in 1:(x$k+2))
      fn.ci_sum[[i]] <- function(y.offset,...) {forestplot::fpDrawSummaryCI(y.offset=0.5,...)}
    if (prediction)
      fn.ci_sum[[x$k+3]] <- function(y.offset,...) {forestplot::fpDrawBarCI(y.offset=0.5,...)}
  }
  # specify colors:
  if (missing(col)) {
    if (shrinkage) {
      col <- forestplot::fpColors(box=c("black", "grey45"),
                                  lines=c("black","grey45"),
                                  summary="grey30")
    } else {
      col <- forestplot::fpColors(box="black", lines="black", summary="grey30")
    }
  }
  # specify plotting sizes:
  if (missing(boxsize)) {
    boxsize <- c(rep(0.25,x$k+1), 0.4)
    if (prediction) boxsize <- c(boxsize, 0.2)
  }
  if (shrinkage && missing(legend))
    legend  <- c("quoted estimate", "shrinkage estimate")
  # specify data for plotting:
  if (shrinkage) { # (show shrinkage intervals)
    mean.arg  <- cbind(ma.dat[,1], ma.shrink[,1])
    lower.arg <- cbind(ma.dat[,2], ma.shrink[,2])
    upper.arg <- cbind(ma.dat[,3], ma.shrink[,3])
  } else {         # (no shrinkage intervals)
    mean.arg  <- ma.dat[,1]
    lower.arg <- ma.dat[,2]
    upper.arg <- ma.dat[,3]
  }
  fp <- NULL
  if (plot) {
    fp <- forestplot::forestplot(labeltext  = labeltext,
                                 mean       = mean.arg,
                                 lower      = lower.arg,
                                 upper      = upper.arg,
                                 is.summary = c(TRUE, rep(FALSE, x$k), TRUE, TRUE),
                                 hrzl_lines = horizl,
                                 fn.ci_norm = fn.ci_norm,
                                 fn.ci_sum  = fn.ci_sum,
                                 col        = col,
                                 boxsize    = boxsize,
                                 legend     = legend, ...)
    plot(fp)
    # add heterogeneity phrase at bottom left:
    if (heterogeneity) {
      tauFigures <- x$summary[c("median","95% lower", "95% upper"), "tau"]
      tauFigures <- tauFigures[tauFigures > 0]
      formatstring <- paste0("%.", decplaces(tauFigures, digits), "f")
      tauphrase <- sprintf(paste0("Heterogeneity (tau): ",formatstring,
                                  " [",formatstring,", ",formatstring,"]"),
                           x$summary["median","tau"],
                           x$summary["95% lower","tau"], x$summary["95% upper","tau"])
      grid::seekViewport("axis_margin")
      tvp <- grid::viewport(x=grid::unit(0.0, "npc"), y=grid::unit(0.0, "npc"),
                            width=grid::stringWidth(tauphrase),
                            height=grid::unit(2, "lines"),
                            just=c("left","bottom"), name="heterogeneityEstimate")
      grid::pushViewport(tvp)
      grid::grid.text(tauphrase, x=grid::unit(0.0, "npc"), y=grid::unit(0.5,"npc"),
                      just=c("left", "centre"), gp=grid::gpar(fontface="oblique"))
    }
  }
  invisible(list("data"       = ma.dat[-1,],
                 "shrinkage"  = ma.shrink[2:(x$k+1),],
                 "labeltext"  = labeltext,
                 "forestplot" = fp))
}



forestplot.escalc <- function(x, labeltext,
                              exponentiate  = FALSE,
                              digits        = 2,
                              plot          = TRUE,
                              fn.ci_norm, fn.ci_sum, col, legend, boxsize, ...)
#
# ARGUMENTS:
#   x            :  a "bayesmeta" object.
#   labeltext    :  you may provide an alternative "labeltext" argument here
#                   (see also the "forestplot()" help).
#   exponentiate :  flag indicating whether to exponentiate numbers (figure and table).
#   prediction   :  flag indicating whether to show prediction interval.
#   shrinkage    :  flag indicating whether to show shrinkage estimates.
#   digits       :  number of significant digits to be shown (based on standard deviations).
#   plot         :  flag you can use to suppress actual plotting.
#   ...          :  further arguments passed to the "forestplot" function.
#   
# VALUE:
#   a list with components
#     $data       :  the meta-analyzed data, and mean and prediction estimates
#     $labeltext  :  the "forestplot()" function's "labeltext" argument used internally
#     $forestplot :  the "forestplot()" function's returned value
#
{
  if (!requireNamespace("forestplot", quietly=TRUE))
    stop("required 'forestplot' package not available!")
  if (utils::packageVersion("forestplot") < "1.5.2")
    warning("you may need to update 'forestplot' to a more recent version (>=1.5.2).")
  # auxiliary function:
  decplaces <- function(x, signifdigits=3)
  # number of decimal places (after decimal point)
  # to be displayed if you want at least "signifdigits" digits
  {
    return(max(c(0, -(floor(log10(x))-(signifdigits-1)))))
  }
  # some sanity checks for the provided arguments:
  stopifnot(is.element("escalc", class(x)),
            length(digits)==1, digits==round(digits), digits>=0,
            length(exponentiate)==1, is.logical(exponentiate),
            length(plot)==1, is.logical(plot))
  
  # extract relevant column names:
  attri <- attributes(x)
  if (all(is.element(c("yi.names", "vi.names"), names(attri)))) { # (for recent "metafor" versions)
    var.names <- c(attri$yi.names, attri$vi.names)
  } else if (is.element("var.names", names(attri))) {             # (for older "metafor" versions)
    var.names <- attri$var.names
  } else {
    stop(paste("Cannont extract \"yi.names\" and \"vi.names\" (or \"var.names\") attribute(s) from escalc object ",
               "(check use of the \"var.names\" option).", sep=""))
  }
  stopifnot(length(var.names)==2, all(is.character(var.names)))
  if (!all(is.element(var.names, names(x)))) {
    stop(paste("Cannont find columns \"",
               var.names[1],"\" and/or \"",
               var.names[2],"\" in escalc object ",
               "(check use of the \"var.names\" option).", sep=""))
  }
  if (is.element("slab", names(attributes(x[,var.names[1]])))) {
    labels <- as.character(attr(x[,var.names[1]], "slab"))
  } else {
    labels <- sprintf("%02d", 1:nrow(x))
  }
  
  # extract data to be plotted:
  y      <- as.vector(x[,var.names[1]])
  sigma  <- sqrt(as.vector(x[,var.names[2]]))
  k <- length(y)
  q95 <- qnorm(0.975)
  ma.dat <- rbind(NA,
                  cbind(y, y - q95*sigma, y + q95*sigma))
  colnames(ma.dat) <- c("estimate", "lower", "upper")
  rownames(ma.dat) <- c("", labels)
  if (exponentiate) {
    ma.dat    <- exp(ma.dat)
  }
  
  # generate "labeltext" data table for plot (unless already provided):
  if (missing(labeltext)) {
    # determine numbers of digits based on standard deviations:
    if (exponentiate) {
      stdevs <- exp(y)*sigma
    } else {
      stdevs <- sigma
    }
    stdevs <- abs(stdevs[is.finite(stdevs) & (stdevs != 0)])
    formatstring <- paste0("%.", decplaces(stdevs, digits), "f")
    # fill data table:
    labeltext <- matrix(NA_character_, nrow=nrow(ma.dat), ncol=3)
    labeltext[1,] <- c("study","estimate", "95% CI")
    labeltext[,1] <- c("study", labels)[1:nrow(ma.dat)]
    for (i in 2:(nrow(ma.dat))) {
      labeltext[i,2] <- sprintf(formatstring, ma.dat[i,"estimate"])
      labeltext[i,3] <- paste0("[", sprintf(formatstring, ma.dat[i,"lower"]),
                                 ", ", sprintf(formatstring, ma.dat[i,"upper"]), "]")
    }
  }
  # add horizontal lines to plot:
  horizl <- list(grid::gpar(col="grey"))
  names(horizl) <- as.character(c(2))
  # specify function(s) for drawing estimates / shrinkage estimates:
  if (missing(fn.ci_norm)) {
    fn.ci_norm <- function(...) {forestplot::fpDrawPointCI(pch=15,...)}
  }
  # specify function(s) for drawing summaries (diamond / bar):
  if (missing(fn.ci_sum)) {
    fn.ci_sum <- list(NULL)
    for (i in 1:(k+1))
      fn.ci_sum[[i]] <- function(y.offset,...) {forestplot::fpDrawSummaryCI(y.offset=0.5,...)}
  }
  # specify colors:
  if (missing(col)) {
    col <- forestplot::fpColors(box="black", lines="black", summary="grey30")
  }
  # specify plotting sizes:
  if (missing(boxsize)) {
    boxsize <- rep(0.25,k+1)
  }
  # specify data for plotting:
  mean.arg  <- ma.dat[,1]
  lower.arg <- ma.dat[,2]
  upper.arg <- ma.dat[,3]
  fp <- NULL
  if (plot) {
    fp <- forestplot::forestplot(labeltext  = labeltext,
                                 mean       = mean.arg,
                                 lower      = lower.arg,
                                 upper      = upper.arg,
                                 is.summary = c(TRUE, rep(FALSE, k)),
                                 hrzl_lines = horizl,
                                 fn.ci_norm = fn.ci_norm,
                                 fn.ci_sum  = fn.ci_sum,
                                 col        = col,
                                 boxsize    = boxsize,
                                 legend     = legend, ...)
    plot(fp)
  }
  invisible(list("data"       = ma.dat[-1,],
                 "labeltext"  = labeltext,
                 "forestplot" = fp))
}



convolve <- function(dens1, dens2,
                     cdf1=Vectorize(function(x){integrate(dens1,-Inf,x)$value}),
                     cdf2=Vectorize(function(x){integrate(dens2,-Inf,x)$value}),
                     delta=0.01, epsilon=0.0001)
# Convolution function from:
#   C. Roever, T. Friede.
#   Discrete approximation of a mixture distribution via restricted divergence.
#   Journal of Computational and Graphical Statistics, 26(1):217-222, 2017.
#   http://doi.org/10.1080/10618600.2016.1276840
#
# (a grid is constructed in the first density's domain ("dens1")
#  so that the convolution is a sum of "dens2" conditionals,
#  weighted by "dens1".)
{
  # some basic sanity checks:
  stopifnot(is.function(dens1), is.function(dens2),
            is.function(cdf1), is.function(cdf2),
            #dens1(-1)>0, dens2(-1)>0,
            delta>0, epsilon>0)
  
  symKL <- function(d)
  # divergence w.r.t. shifting of the 2nd distribution ("dens2()")
  {
    stopifnot(d>=0)
    func1 <- function(x)
    # integrand for one directed divergence
    {
      d2md <- dens2(x-d)
      # utilize density's "log" argument, if possible:
      if (is.element("log", names(formals(dens2)))) {
        logd2   <- dens2(x, log=TRUE)
        logd2md <- dens2(x-d, log=TRUE)
      } else {
        logd2   <- log(dens2(x))
        logd2md <- log(d2md)
      }
      return(ifelse((logd2>-Inf) & (logd2md>-Inf), (logd2md-logd2)*d2md, 0.0))
    }
    
    func2 <- function(x)
    # integrand for other directed divergence
    {
      d2 <- dens2(x)
      # utilize density's "log" argument, if possible:
      if (is.element("log", names(formals(dens2)))) {
        logd2   <- dens2(x, log=TRUE)
        logd2md <- dens2(x-d, log=TRUE)
      } else {
        logd2   <- log(d2)
        logd2md <- log(dens2(x-d))
      }
      return(ifelse((logd2>-Inf) & (logd2md>-Inf), (logd2-logd2md)*d2, 0.0))
    }
    int1 <- integrate(func1, -Inf, Inf)
    if (int1$message != "OK")
      warning(paste0("Problem computing KL-divergence (1): \"", int1$message,"\""))
    int2 <- integrate(func2, -Inf, Inf)
    if (int2$message != "OK")
      warning(paste0("Problem computing KL-divergence (2): \"", int2$message,"\""))
    return(int1$value + int2$value)
  }

  # determine bin half-width:
  step <- sqrt(delta)
  while (symKL(step) < delta) step <- 2*step
  ur <- uniroot(function(X){return(symKL(X)-delta)}, lower=0, upper=step)
  step <- ur$root
  # determine grid range via  epsilon/2  and  1-epsilon/2  quantiles:
  mini <- -1
  while (cdf1(mini) > epsilon/2) mini <- 2*mini
  maxi <- 1
  while (cdf1(maxi) < 1-(epsilon/2)) maxi <- 2*maxi
  ur <- uniroot(function(X){return(cdf1(X)-epsilon/2)},
                lower=mini, upper=maxi)
  mini <- ur$root
  ur <- uniroot(function(X){return(cdf1(X)-(1-epsilon/2))},
                lower=mini, upper=maxi)
  maxi <- ur$root
  # determine number of reference points:
  k <- ceiling((maxi-mini)/(2*step))+1
  # determine reference points:
  support <- mini - ((k*2*step)-(maxi-mini))/2 + (0:(k-1))*2*step
  # determine bin margins:
  margins <- support[-1]-step
  # determine bin weights:
  weight <- rep(NA, length(support))
  for (i in 1:(k-1))
    weight[i] <- cdf1(margins[i])
  weight[k] <- 1
  for (i in k:2)
    weight[i] <- weight[i]-weight[i-1]
  # the eventual grid:
  grid <- cbind("lower"=c(-Inf, margins),
                "upper"=c(margins, Inf),
                "reference"=support,
                "prob"=weight)

  # probability density of convolution:
  density <- function(x)
  {
    return(apply(matrix(x,ncol=1), 1,
                 function(x){sum(grid[,"prob"]*dens2(x-grid[,"reference"]))}))
  }

  # cumulative distribution function (CDF) of convolution:
  cdf <- function(x)
  {
    return(apply(matrix(x,ncol=1), 1,
                 function(x){sum(grid[,"prob"]*cdf2(x-grid[,"reference"]))}))
  }

  # quantile function (inverse CDF) of convolution:
  quantile <- function(p)
  {
    quant <- function(pp)
    {
      mini <- -1
      while (cdf(mini) > pp) mini <- 2*mini
      maxi <- 1
      while (cdf(maxi) < pp) maxi <- 2*maxi
      ur <- uniroot(function(x){return(cdf(x)-pp)}, lower=mini, upper=maxi)
      return(ur$root)      
    }
    proper <- ((p>0) & (p<1))
    result <- rep(NA,length(p))
    if (any(proper)) result[proper] <- apply(matrix(p[proper],ncol=1), 1, quant)
    return(result)
  }
  
  return(list("delta"    = delta,     # tuning parameter (maximum divergence)
              "epsilon"  = epsilon,   # tuning parameter (neglected tail probability)
              "binwidth" = 2*step,    # bin width (constant across bins)
              "bins"     = k,         # number of bins
              "support"  = grid,      # bin margins, reference points and weights
              "density"  = density,   # resulting probability density function
              "cdf"      = cdf,       # resulting cumulative distribution function (CDF)
              "quantile" = quantile)) # resulting quantile function (inverse CDF)
}



funnel.bayesmeta <- function(x,
                             main=deparse(substitute(x)),
                             xlab=expression("effect "*y[i]),
                             ylab=expression("standard error "*sigma[i]),
                             zero=0.0, FE=FALSE, legend=FE, ...)
# generate funnel plot from a "bayesmeta" object.
{
  # sanity check:
  stopifnot(is.element("bayesmeta", class(x)))
  # range of standard error axis:
  yrange <- c(0.0, max(x$sigma))
  # standard error levels for prediction intervals:
  sevec <- seq(from=0, yrange[2]*1.04, le=27)
  # compute (RE) prediction intervals
  # (for observed "y" values, conditional on standard error):
  intRE <- matrix(NA_real_, nrow=length(sevec), ncol=2,
                dimnames=list(NULL, c("lower","upper")))
  intRE[1,] <- x$qposterior(theta.p=c(0.025, 0.975), predict=TRUE)
  for (i in 2:length(sevec)){
    conv <- try(convolve(dens1=function(a, log=FALSE){return(x$dposterior(theta=a, predict=TRUE, log=log))},
                         dens2=function(b, log=FALSE){return(dnorm(x=b, mean=0, sd=sevec[i], log=log))},
                         cdf1 =function(a){return(x$pposterior(theta=a, predict=TRUE))},
                         cdf2 =function(b){return(pnorm(q=b, mean=0, sd=sevec[i]))}))
    if (all(class(conv)!="try-error")) {
      intRE[i,] <- conv$quantile(p=c(0.025, 0.975))
    }
  }
  # compute _FE_ prediction intervals
  # (for observed "y" values, conditional on standard error):
  intFE <- matrix(NA_real_, nrow=length(sevec), ncol=2,
                  dimnames=list(NULL, c("lower","upper")))
  cm <- x$cond.moment(tau=0)
  for (i in 1:length(sevec)){
    intFE[i,] <- qnorm(c(0.025, 0.975), mean=cm[1,"mean"], sd=sqrt(cm[1,"sd"]^2+sevec[i]^2))
  }
  FEcol="red3"
  REcol="blue3"
  # empty plot:
  plot(range(intRE), -yrange, type="n",
       ylab=ylab, xlab=xlab, main=main, axes=FALSE)
  # funnels (grey area):
  polygon(c(intRE[,1], rev(intRE[,2])), c(-sevec, rev(-sevec)), col="grey90", border=NA)
  if (FE) polygon(c(intFE[,1], rev(intFE[,2])), c(-sevec, rev(-sevec)), col="grey80", border=NA)
  # grid lines:
  lines(c(intRE[1,1], intRE[1,1], NA, intRE[1,2], intRE[1,2]),
        c(0,-max(sevec), NA, 0, -max(sevec)), col="grey75", lty="dashed")
  abline(h=0, col="darkgrey")
  yticks <- pretty(yrange)
  abline(h=-yticks[yticks>0], col="grey75", lty="15")
  # funnels (outline):
  matlines(intRE, cbind(-sevec, -sevec), col=REcol, lty="dashed")
  if (FE) matlines(intFE, cbind(-sevec, -sevec), col=FEcol, lty="dotted")
  lines(rep(x$summary["median","theta"], 2), range(-sevec), col=REcol, lty="dashed")
  if (FE) lines(rep(cm[1,"mean"], 2), range(-sevec), col=FEcol, lty="dotted")
  # zero line:
  if (is.finite(zero))
    lines(c(zero, zero), c(-1,1)*max(sevec), col="darkgrey", lty="solid")
  # actual points:
  points(x$y, -x$sigma, pch=21, col="black", bg=grey(0.5, alpha=0.5), cex=1)
  if (FE && legend)
    legend("topleft", c("RE model", "FE model"),
           col=c(REcol, FEcol), lty=c("dashed", "dotted"), bg="white")
  axis(1); axis(2, at=-yticks, labels=yticks); box()
  invisible()
}



normalmixture <- function(density,
                          cdf = Vectorize(function(x){integrate(density,0,x)$value}),
                          mu = 0,
                          delta = 0.01, epsilon = 0.0001,
                          rel.tol.integrate=2^16*.Machine$double.eps,
                          abs.tol.integrate=rel.tol.integrate,
                          tol.uniroot=rel.tol.integrate)
# derive normal mixture where mean (mu) is given (and fixed)
# and the standard deviation (sigma) follows a distribution
# given through "density" or "cdf" argument.
# Mixture approximation is done using the 'DIRECT' algorithm
# described in http://arxiv.org/abs/1602.04060
{
  # some basic sanity checks:
  stopifnot((!missing(cdf) || !missing(density)),
            is.function(cdf), delta>0, epsilon>0)
  if (missing(cdf) && !missing(density)) {
    stopifnot(is.function(density))
    # check for properness of mixing distribution:
    density.integral <- integrate(density, lower=0, upper=Inf,
                                  rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)
    stopifnot(density.integral$message == "OK",
              (abs(density.integral$value - 1.0) <= density.integral$abs.error))
  }

  nextsigma <- function(sigma=1, delta=0.01)
  # analytical solution for next larger sigma value
  # for which KL-divergence w.r.t. current sigma is = delta
  {
     return(sigma * sqrt(1 + delta^2 + sqrt(delta^2 + 2*delta)))
  }

  # determine (minimally required) grid range:
  s <- 1
  while (cdf(s) <= epsilon/2) s <- s*2
  ur <- uniroot(function(x){cdf(x)-epsilon/2}, lower=0, upper=s, tol=tol.uniroot)
  lower <- ur$root
  while (cdf(s) <= 1-epsilon/2) s <- s*2
  ur <- uniroot(function(x){cdf(x)-(1-epsilon/2)}, lower=lower, upper=s, tol=tol.uniroot)
  upper <- ur$root

  # now fill grid:
  grid <- matrix(c(0, NA, lower, NA),
                 nrow=1,
                 dimnames=list(NULL, c("lower","upper","reference","prob")))
  s <- lower
  i <- 1
  while (s < upper) {
    i <- i+1
    grid <- rbind(grid, NA)
    s <- nextsigma(s, delta=delta)
    grid[i-1,"upper"] <- grid[i,"lower"] <- s
    s <- nextsigma(s, delta=delta)
    grid[i,"reference"] <- s
  }
  grid[i, "upper"] <- Inf

  grid[,"prob"] <- diff(c(0, cdf(grid[,"upper"])))
  grid[,"prob"] <- grid[,"prob"] / sum(grid[,"prob"])

  # probability density of mixture:
  dens <- function(x)
  {
    return(apply(matrix(x,ncol=1), 1,
                 function(y){sum(grid[,"prob"]*dnorm(y, mean=mu, sd=grid[,"reference"]))}))
  }

  # cumulative distribution function (CDF) of mixture:
  cumul <- function(x)
  {
    return(apply(matrix(x,ncol=1), 1,
                 function(y){sum(grid[,"prob"]*pnorm(y, mean=mu, sd=grid[,"reference"]))}))
  }

  # quantile function (inverse CDF) of mixture:
  quantile <- function(p)
  {
    quant <- function(pp)
    {
      mini <- mu-1
      while (cumul(mini) > pp) mini <- mu - 2*(mu-mini)
      maxi <- mu+1
      while (cumul(maxi) < pp) maxi <- mu + 2*(maxi-mu)
      ur <- uniroot(function(x){return(cumul(x)-pp)}, lower=mini, upper=maxi, tol=tol.uniroot)
      return(ur$root)      
    }
    proper <- ((p>0) & (p<1))
    result <- rep(NA,length(p))
    if (any(proper)) result[proper] <- apply(matrix(p[proper],ncol=1), 1, quant)
    return(result)
  }  
  
  result <- list("delta"          = delta,
                 "epsilon"        = epsilon,
                 "mu"             = mu,
                 "bins"           = i,
                 "support"        = grid,
                 "density"        = dens,      # resulting probability density function
                 "cdf"            = cumul,     # resulting cumulative distribution function (CDF)
                 "quantile"       = quantile,  # resulting quantile function (inverse CDF)
                 "mixing.density" = NULL,
                 "mixing.cdf"     = cdf)  
  if (!missing(density)) result$mixing.density <- density
  return(result)
}


pppvalue<- function(x,
                    parameter = "mu",
                    value = 0.0,
                    alternative = c("two.sided", "less", "greater"),
                    statistic = "median",  # a.k.a. "discrepancy variable"
                    rejection.region,
                    n = 10,
                    prior = FALSE,
                    quietly = FALSE,
                    parallel, seed, ...)
# Posterior predictive p-values
#    * Gelman & al., Bayesian Data Analysis, 3rd edition, Chapter 6,
#      Chapman & Hall / CRC, Boca Raton, 2014.
#    * Meng, X.-L. Posterior predictive p-values.
#      The Annals of Statistics, 22(3):1142-1160, 1994.
#      http://doi.org/10.1214/aos/1176325622
#
# Parameters:
#   x                :  a "bayesmeta" object
#   parameter        :  the parameter to be tested
#   value            :  the null-hypothesized value
#   alternative      :  the alternative to be tested against
#   statistic        :  the figure to be used as "test statistic"
#   rejection.region :  the test statistic's rejection region. May be
#                       one of "upper.tail", "lower.tail" or "two.tailed".
#                       If unspecified, it is set based on the "alternative" parameter
#   n                :  the number of Monte Carlo samples
#   prior            :  flag to request _PRIOR_ predictive p-values
#   quietly          :  flag to indicate command line text output
#   parallel         :  number of parallel processes to use
#   ...              :  further arguments passed to "statistic",
#                       if "statistic" argument is a function
#
{
  # some preliminary sanity checks;
  # "x" argument:
  if (!is.element("bayesmeta",class(x)))
    warning("function applicable to objects of class 'bayesmeta' only.")
  stopifnot(is.element("bayesmeta",class(x)))
  # "parameter" argument:
  stopifnot(length(parameter)==1)
  if (is.numeric(parameter)){
    stopifnot(is.element(parameter, 1:x$k))
    parameter <- x$labels[parameter]
  }
  parameter <- match.arg(parameter, c("mu", "tau", x$labels))
  thetapar <- FALSE
  if (!is.element(parameter, c("mu","tau"))) {
    thetapar <- TRUE  # indicates that parameter concerns one of the "theta" (shrinkage) parameters
    indiv.which <- which(is.element(x$labels, parameter))  # index of concerned "theta" parameter
  }
  # "value" argument:
  stopifnot(length(value)==1, is.finite(value),
            (parameter != "tau") || (value >= 0.0))
  # "alternative" argument:
  alternative <- match.arg(alternative)
  # "statistic" argument:
  statfun <- is.function(statistic)
  statNA <- (!statfun && is.na(statistic))
  stopifnot(statfun | is.character(statistic) | statNA)
  if (statNA) {
    statname <- "statistic"
  } else if (statfun) {
    statname <- deparse(substitute(statistic))
  } else {
    statistic <- match.arg(tolower(statistic), c("t", "q", "cdf", rownames(x$summary)))
    if ((parameter=="tau") && (value==0.0)) {
      if (alternative != "greater")
        warning(paste0("a combination of 'parameter=\"tau\"' ",
                       "and 'alternative=\"", alternative, "\"' may not make sense!"))
      if (statistic == "cdf")
        warning(paste0("a combination of 'parameter=\"tau\"', ",
                       "'value=0.0' and 'statistic=\"cdf\"' does not make sense!"))
    }
    # try to speed up 'bayesmeta()' computations by omitting unnecessary CI optimizations:
    inttype <- ifelse(is.element(statistic, c("95% lower", "95% upper")),
                      x$interval.type, "central")
    statname <- ifelse(statistic=="q", "Q", statistic)
  }
  # "rejection.region" argument:
  if (!missing(rejection.region)) {
    rejection.region <- match.arg(rejection.region, c("upper.tail", "lower.tail", "two.tailed"))
  } else { # decide on rejection region based on hypothesis:
    if (statNA | (!statfun && is.element(statistic, c("q"))))
      rejection.region <- "upper.tail"
    else if (!statfun && is.element(statistic, c("cdf")))
      rejection.region <- switch(alternative,
                                 two.sided = "two.tailed",
                                 less      = "upper.tail",
                                 greater   = "lower.tail")
    else
      rejection.region <- switch(alternative,
                                 two.sided = "two.tailed",
                                 less      = "lower.tail",
                                 greater   = "upper.tail")
  }
  stopifnot(is.element(rejection.region, c("upper.tail", "lower.tail", "two.tailed")))
  # "prior" argument:
  if (prior && (!x$tau.prior.proper || !all(is.finite(x$mu.prior))))
    warning("prior predictive p-values require proper priors for effect and heterogeneity!")
  if (prior && !is.element(parameter, c("mu", "tau")))
    warning("prior predictive p-values are only available for effect (mu) and heterogeneity (tau) parameters!")
  stopifnot(!prior || (x$tau.prior.proper && all(is.finite(x$mu.prior))))
  # determine a sensible number of parallel processes:
  if (missing(parallel)) {
    # by default, use "parallel" package if available:
    if (requireNamespace("parallel")) {
      # by default, use all but one core:
      parallel <- max(c(1, parallel::detectCores() - 1))
      # but don't use more cores than MCMC samples:
      parallel <- min(c(parallel, n))
    } else {
      parallel <- 1
    }
  } else { # otherwise check number of cores & processes:
    stopifnot(parallel >= 1)
    if (requireNamespace("parallel")) {
      if (parallel > parallel::detectCores())
        warning("number of requested parallel processes (", parallel,
                ") is larger than number of cores (", parallel::detectCores(), ").")
      if (parallel > n)
        warning("number of requested parallel processes (", parallel,
                ") is larger than number of MC samples (", n, ").")
    } else {
      warning("failed to load \"parallel\" package.")
    }
  }
  seed.missing <- missing(seed)
  stopifnot(seed.missing || (seed==round(seed)))
  ptm <- proc.time()
  ##################
  if (prior) {         # ...conditional prior p(tau | mu)
    ptau <- function(tau)
    # cumulative distribution function (CDF) of tau prior
    {
      stopifnot(length(tau)==1, is.finite(tau), tau>=0)
      return(integrate(function(t){x$dprior(tau=t)}, lower=0, upper=tau,
                       rel.tol=x$rel.tol.integrate, abs.tol=x$abs.tol.integrate)$value)
    }
    ptau <- Vectorize(ptau)
    qtau <- function(p.tau)
    # quantile function (inverse CDF) of tau prior
    {
      stopifnot(length(p.tau)==1, is.finite(p.tau), p.tau>=0, p.tau<=1)
      if (p.tau==0) result <- 0.0
      else if (p.tau==1) result <- Inf
      else {
        upper <- 1.0
        while (ptau(upper) < p.tau) upper <- 2*upper
        result <- uniroot(function(t){ptau(tau=t)-p.tau},
                          lower=0, upper=upper, tol=x$tol.uniroot)$root
      }
      return(result)
    }
  }
  if (parameter=="mu") { # define functions to generate draws from conditional (tau|mu):
    # lookup table for tau conditionals' normalizing constants:
    normconst <- matrix(nrow=0, ncol=2, dimnames=list(NULL, c("mu","const")))
    dctau <- function(tau, mu)
    # conditional density of (tau | mu)
    {
      stopifnot(length(tau)==1, is.finite(tau), tau>=0)
      if ((nrow(normconst)>0) && is.element(mu, normconst[,"mu"])) {
        const <- normconst[which(normconst[,"mu"]==mu),"const"]
      } else {
        const <- integrate(function(t){x$dposterior(tau=t, mu=mu)},
                           lower=0, upper=Inf,
                           rel.tol=x$rel.tol.integrate, abs.tol=x$abs.tol.integrate)$value
        normconst <<- rbind(normconst, c(mu,const)) # (change matrix GLOBALLY)
      }
      return(x$dposterior(tau=tau, mu=mu) / const)
    }
    dctau <- Vectorize(dctau)
    pctau <- function(tau, mu)
    # conditional cumulative distribution function (CDF) of (tau | mu)
    {
      stopifnot(length(tau)==1, is.finite(tau), tau>=0)
      return(integrate(function(t){dctau(tau=t, mu=mu)}, lower=0, upper=tau,
                       rel.tol=x$rel.tol.integrate, abs.tol=x$abs.tol.integrate)$value)
    }
    pctau <- Vectorize(pctau)
    qctau <- function(p.tau, mu)
    # conditional quantile function (inverse CDF) of (tau | mu)
    {
      stopifnot(length(p.tau)==1, is.finite(p.tau), p.tau>=0, p.tau<=1)
      if (p.tau==0) result <- 0.0
      else if (p.tau==1) result <- Inf
      else {
        upper <- 1.0
        while (pctau(upper,mu) < p.tau) upper <- 2*upper
        result <- uniroot(function(t){pctau(tau=t,mu=mu)-p.tau},
                          lower=0, upper=upper, tol=x$tol.uniroot)$root
      }
      return(result)
    }
    rctau <- function(mu)
    # random number generation for (tau | mu)
    {
      return(qctau(p.tau=runif(n=1), mu=mu))
    }
  } else if (thetapar) {
    # compute conditional posterior, conditional on theta[i]==value:
    condy     <- x$y
    condsigma <- x$sigma
    condy[indiv.which]     <- value
    condsigma[indiv.which] <- 0.0
    if (alternative == "two.sided") {
      condbm <- bayesmeta(y=condy, sigma=condsigma, labels=x$labels,
                          mu.prior=x$mu.prior,
                          tau.prior=x$dprior,
                          interval.type=inttype,
                          rel.tol.integrate=x$rel.tol.integrate,
                          abs.tol.integrate=x$abs.tol.integrate,
                          tol.uniroot=x$tol.uniroot)
    }
  }
  # determine realized "test statistic" value in actual data:
  if (statNA) {
    stat.actual <- NA_real_
  } else if (statfun) {
    stat.actual <- statistic(x$y, ...)
  } else {
    if (statistic == "q") {
      cochranQ <- function(yi, si)
      {
        wi <- 1/si^2
        yhat <- sum(yi * wi) / sum(wi)
        Q <- sum(((yi-yhat)/si)^2)
        return(Q)
      }
      stat.actual <- cochranQ(x$y, x$sigma)
    } else if (statistic=="cdf") {
      if (parameter=="tau")
        stat.actual <- x$pposterior(tau=value)
      else if (parameter=="mu")
        stat.actual <- x$pposterior(mu=value)
      else
        stat.actual <- x$pposterior(theta=value, individual=indiv.which)
    } else {
      if (thetapar) {
        if (statistic == "t")
          stat.actual <- (x$theta["mean", indiv.which] - value) / x$theta["sd", indiv.which]
        else
          stat.actual <- x$theta[statistic, indiv.which]
      } else {
        if (statistic == "t")
          stat.actual <- (x$summary["mean", parameter] - value) / x$summary["sd", parameter]
        else
          stat.actual <- x$summary[statistic, parameter]
      }
    }
  }
  names(stat.actual) <- statname
  # determine at which (prior/posterior) quantile hypothesized value is situated
  # (necessary for sampling for single-tailed hypotheses):
  if (alternative != "two.sided") {
    if (parameter == "mu") {
      if (prior) p.mu <- pnorm(q=value, mean=x$mu.prior[1], sd=x$mu.prior[2])  # prior quantile
      else       p.mu <- x$pposterior(mu=value)                                # posterior quantile
    } else if (parameter == "tau") {
      if (prior) p.tau <- ptau(value)              # prior quantile
      else       p.tau <- x$pposterior(tau=value)  # posterior quantile
    } else {
      p.theta <- x$pposterior(theta=value, individual=indiv.which)  # posterior quantile
    }
  }
  
  # the main function:
  pppfun <- function(i)
  {
    if (! seed.missing) set.seed(seed + i)
    # generate data:
    sigma <- x$sigma
    if (parameter=="mu") {          # Null hypothesis concerns effect mu
      if (alternative=="two.sided") {   # fix effect mu at hypothesized value:
        rmu <- value
      } else if (alternative=="less") { # draw effect mu from H0's domain's conditional:
        if (prior) rmu <- qnorm(runif(1, p.mu, 1.0), mean=x$mu.prior[1], sd=x$mu.prior[2])  # prior
        else       rmu <- x$qposterior(mu.p=runif(1, p.mu, 1.0))                            # posterior
      } else {
        if (prior) rmu <- qnorm(runif(1, 0.0, p.mu), mean=x$mu.prior[1], sd=x$mu.prior[2])  # prior
        else       rmu <- x$qposterior(mu.p=runif(1, 0.0, p.mu))                            # posterior
      }
      # draw tau from conditional (tau|mu):
      if (prior) rtau <- qtau(runif(1, 0.0, 1.0))  # prior
      else       rtau <- rctau(mu=rmu)             # posterior
      # draw study-specific effects (theta) and effects (y):
      rtheta <- rnorm(n=x$k, mean=rmu,    sd=rtau)
      y      <- rnorm(n=x$k, mean=rtheta, sd=sigma)
    } else if (parameter=="tau") {  # Null hypothesis concerns heterogeneity tau
      if (alternative=="two.sided") {   # fix heterogeneity tau at hypothesized value:
        rtau <- value
      } else if (alternative=="less") { # draw effect mu from H0's domain's conditional:
        if (prior) rtau <- qtau(runif(1, p.tau, 1.0))
        else       rtau <- x$qposterior(tau.p=runif(1, p.tau, 1.0))
      } else {
        if (prior) rtau <- qtau(runif(1, 0.0, p.tau))
        else       rtau <- x$qposterior(tau.p=runif(1, 0.0, p.tau))
      }
      # draw mu from conditional (mu|tau):
      if (prior) {
        rmu <- qnorm(runif(1, 0.0, 1.0), mean=x$mu.prior[1], sd=x$mu.prior[2])
      } else {
        cm <- x$cond.moment(tau=rtau)
        rmu <- rnorm(n=1, mean=cm[1,"mean"], sd=cm[1,"sd"])
      }
      # draw study-specific effects (theta) and effects (y):
      rtheta <- rnorm(n=x$k, mean=rmu,    sd=rtau)
      y      <- rnorm(n=x$k, mean=rtheta, sd=sigma)
    } else {                        # Null hypothesis concerns ith "shrinkage" parameter theta
      if (alternative=="two.sided"){
        rtheta.i <- value
        taumu <- condbm$rposterior(n=1)[1,]
      } else {
        if (alternative=="less") {
          rtheta.i <- x$qposterior(theta.p=runif(1, p.theta, 1.0), indiv=indiv.which)
        } else {
          rtheta.i <- x$qposterior(theta.p=runif(1, 0.0, p.theta), indiv=indiv.which)
        }
        condy[indiv.which] <- rtheta.i
        condbm <- bayesmeta(y=condy, sigma=condsigma, labels=x$labels,
                            mu.prior=x$mu.prior,
                            tau.prior=x$dprior,
                            interval.type=inttype,
                            rel.tol.integrate=x$rel.tol.integrate,
                            abs.tol.integrate=x$abs.tol.integrate,
                            tol.uniroot=x$tol.uniroot)
        taumu <- condbm$rposterior(n=1)[1,]
      }
      rtau <- taumu["tau"]
      rmu <- taumu["mu"]
      rtheta <- rnorm(n=x$k, mean=rmu, sd=rtau)
      rtheta[indiv.which] <- rtheta.i
      #y      <- rnorm(n=x$k, mean=rmu, sd=sigma)
      y      <- rnorm(n=x$k, mean=rtheta, sd=sigma)
    }
    # data (y) generated. Now compute statistic:
    if (statNA) {
      statval <- NA_real_
    } else if (statfun) {
      statval <- statistic(y, ...)
    } else {
      if (statistic == "q") {
        statval <- cochranQ(y, x$sigma)
      } else {
        # perform meta-analysis:
        bm <- bayesmeta(y=y, sigma=x$sigma, labels=x$labels,
                        mu.prior=x$mu.prior,
                        tau.prior=x$dprior,
                        interval.type=inttype,
                        rel.tol.integrate=x$rel.tol.integrate,
                        abs.tol.integrate=x$abs.tol.integrate,
                        tol.uniroot=x$tol.uniroot)
        if (thetapar) {
          if (statistic == "t")
            statval <- (bm$theta["mean", indiv.which] - value) / bm$theta["sd", indiv.which]
          else if (statistic=="cdf")
            statval <- bm$pposterior(theta=value, individual=indiv.which)
          else
            statval <- bm$theta[statistic, indiv.which]
        } else {
          if (statistic == "t")
            statval <- (bm$summary["mean", parameter] - value) / bm$summary["sd", parameter]
          else if (statistic=="cdf") {
            if (parameter=="tau")
              statval <- bm$pposterior(tau=value)
            else if (parameter=="mu")
              statval <- bm$pposterior(mu=value)
          } else
            statval <- bm$summary[statistic, parameter]
        }
      }
    }
    # statistic computed. Return results:
    result <- unname(c(rtau, rmu, statval, rtheta, y))
    names(result) <- c("tau", "mu", "statistic",
                       sprintf("theta[%d]", 1:x$k),
                       sprintf("y[%d]", 1:x$k))
    return(result)
  } # END pppfun()
    
  # break calculations down into smaller bits
  # to allow for progress bar updates (unless suppressed)
  if (quietly) # (do computations in one go)
    idxlist <- list(1:n)
  else {
    # determine intermediate breakpoints (to update progress bar)
    # at multiples of 'parallel':
    stopfrac <- c(0.01, 0.02, 0.05, (1:9)/10)  # (1%, 2%, 5%, 10%, ...)
    stopint <- unique(ceiling(stopfrac * (n/parallel)))
    stopint <- c(stopint * parallel, n)
    stopint <- unique(stopint[stopint <= n])
    # assemble list of indices to process at each stage:
    idxlist <- list(1:stopint[1])
    if (length(stopint) > 1)
      for (i in 2:length(stopint))
        idxlist[[i]] <- (max(idxlist[[i-1]])+1):stopint[i]
  }
  
  #  the "hot loop":
  stat.repli <- NULL
  if (parallel > 1) { # initialize cluster:
    clust <- parallel::makeCluster(parallel)
    #parallel::clusterEvalQ(clust, library("bayesmeta", lib.loc="~/temp/test"))  #  <--  /!\  HACK!
    parallel::clusterEvalQ(clust, library("bayesmeta"))
  }
  if (!quietly) {      # (initialize progress bar etc.)
    cat(paste0("  Generating n=",n," Monte Carlo samples.\n"))
    if (n <= 100)
      cat(paste0("  /!\\  Caution: a sample size of  n >> 100  will usually be appropriate.\n"))
    cat(paste0("  Sampling progress",
               ifelse(parallel > 1, paste0(" (using ",parallel," parallel processes)"), ""),
               ":\n"))
    pb <- utils::txtProgressBar(style=3)
    utils::setTxtProgressBar(pb, 0.0)
  }
  for (i in 1:length(idxlist)){
    if (parallel == 1) { # simple "sapply()" call
      stat.repli <- rbind(stat.repli,
                          t(sapply(idxlist[[i]], pppfun, simplify=TRUE)))
    } else {             # parallel computation via "parSapply()":
      stat.repli <- rbind(stat.repli,
                          t(parallel::parSapply(clust, idxlist[[i]], pppfun, simplify=TRUE)))
    }
    if (!quietly) utils::setTxtProgressBar(pb, max(idxlist[[i]])/n)
  }
  if (parallel > 1) { # terminate cluster:
    parallel::stopCluster(clust)
  }
  
  # how many replicates are in the "tails" beyond the actualized statistic value:
  tail <- factor(rep("no", nrow(stat.repli)),
                 levels=c("lower", "no", "upper"), ordered=TRUE)
  lowertail <- (stat.repli[,"statistic"] <= stat.actual)
  uppertail <- (stat.repli[,"statistic"] >= stat.actual)
  if (rejection.region=="lower.tail")
    tail[lowertail] <- "lower"
  else if (rejection.region=="upper.tail")
    tail[uppertail] <- "upper"
  else {  # note: two-sided version is based on equal-tailed rejection region
    if (sum(uppertail, na.rm=TRUE) < sum(lowertail, na.rm=TRUE)) {
      tail[uppertail] <- "upper"
      tail[rank(stat.repli[,"statistic"], ties.method="min") <= sum(uppertail, na.rm=TRUE)] <- "lower"
    } else {
      tail[lowertail] <- "lower"
      tail[rank(stat.repli[,"statistic"], ties.method="max") > (n-sum(lowertail, na.rm=TRUE))] <- "upper"
    }
  }
  tail[is.na(stat.repli[,"statistic"])] <- NA
  n.tail <- sum(tail != "no", na.rm=TRUE) + sum(is.na(tail))
  # corresponding p-value:
  quant <- ifelse(statNA, NA_real_, n.tail / n)
  # null hypothesis:
  nullval <- value
  names(nullval) <- ifelse(thetapar,
                           paste0("study effect (",indiv.which,": '",x$labels[indiv.which],"')"),
                           ifelse(parameter=="mu", "effect (mu)", "heterogeneity (tau)"))
  replicates <- list("tau"       = stat.repli[,"tau"],
                     "mu"        = stat.repli[,"mu"],
                     "theta"     = stat.repli[,(3+(1:x$k))],
                     "y"         = stat.repli[,(3+x$k+(1:x$k))],
                     "statistic" = cbind.data.frame(stat.repli[,"statistic"], "tail"=tail))
  colnames(replicates$theta) <- colnames(replicates$y) <- x$labels
  colnames(replicates$statistic) <- c(statname, "tail")
  
  ptm <- proc.time() - ptm
  if (!quietly) cat(paste0("\n  (computation time: ",
                          sprintf("%.1f",ptm[3]), " seconds = ",
                          sprintf("%.1f",ptm[3]/60), " minutes.)\n"))
  
  result <- list("statistic"        = stat.actual,
                 "parameter"        = c("Monte Carlo replicates"=n),
                 "p.value"          = quant,
                 "null.value"       = nullval,
                 "alternative"      = alternative,
                 "method"           = paste0("'bayesmeta' ",
                                             ifelse(prior, "prior", "posterior"), " predictive p-value (",
                                             ifelse(alternative=="two.sided", "two-sided", "one-sided"), ")"),
                 "data.name"        = deparse(substitute(x)),
                 # nonstandard-"htest"-elements:
                 "call"             = match.call(expand.dots=FALSE),
                 "rejection.region" = rejection.region,
                 "replicates"       = replicates,
                 "computation.time" = c("seconds"=unname(ptm[3])))
  class(result) <- "htest"
  return(result)
}



uisd <- function(n, ...)
{
  UseMethod("uisd")
}


uisd.default <- function(n, sigma, sigma2=sigma^2, labels=NULL, individual=FALSE, ...)
# Compute unit information standard deviation (UISD).
# Arguments:
#   nc         :  studies' sample sizes
#   sigma      :  studies' standard errors
#   sigma2     :  studies' variances (squared standard errors)
#   individual :  if `TRUE', study-specific UISDs are returned
#                 (instead of overall)
{
  stopifnot(length(n)==length(sigma2),
            all(is.finite(n)), all(is.finite(sigma2)),
            all(n>0), all(sigma2>0),
            (!individual) |
             ((length(labels)==0) || (length(labels)==length(sigma2))))
  
  if (individual) {
    if (!is.character(labels))
      labels <- as.character(labels)
    result <- sqrt(n*sigma2)
    names(result) <- labels
  } else {
    result <- sqrt(sum(n) / sum(1/sigma2))
  }
  return(result)
}


uisd.escalc <- function(n, ...)
# Compute unit information standard deviation (UISD).
# Arguments:
#   n :  an "escalc" object
{
  stopifnot(is.element("escalc", class(n)))
  if (!is.element("ni", names(attributes(n$yi))))
    stop("An \"ni\" attribute is required for the \"escalc\" object!")
  return(uisd(n=attr(n$yi, "ni"), sigma2=n$vi,
              labels=attr(n$yi, "slab"), ...))
}


ess <- function(object, ...)
{
  UseMethod("ess")
}


ess.bayesmeta <- function(object, uisd,
                          method=c("elir", "vr", "pr", "mtm.pt"), ...)
{
  # Computation of effective sample sizes; see:
  #
  #   B. Neuenschwander, S. Weber, H. Schmidli, A. O'Hagan.
  #   Predictively consistent prior effective sample sizes.
  #   Biometrics 76(2): 578-587, 2020.
  #   https://doi.org/10.1111/biom.13252
  #
  # Note that  i_F(theta)  here equals  1 / UISD^2
  # (where the UISD may or may not vary with theta).

  # preliminary sanity checks for provided arguments:
  stopifnot(is.element("bayesmeta", class(object)))
  if (missing(uisd)) {
    stop("\"uisd\" argument need to be specified (either a numeric value or a function)!")
  }
  stopifnot(is.vector(uisd) | is.function(uisd))
  if (is.vector(uisd)) { # specify function returning a constant:
    stopifnot(length(uisd)==1, is.numeric(uisd),
              is.finite(uisd), uisd > 0)
    uisdFun <- function(theta) {return(rep(uisd, length(theta)))}
  } else {               # specify a wrapper function (with built-in sanity checks):
    uisdFun <- function(theta)
    {
      result <- uisd(theta)
      if (!is.vector(result) | !is.numeric(result)) {
        warning(paste0("Invalid output from \"uisd()\" function: uisd(",theta,")"))
      }
      if (!is.finite(result)) {
        warning(paste0("\"uisd()\" function needs to return finite values! (uisd(",theta,") <= 0)"))
      }
      if (result <= 0.0) {
        warning(paste0("\"uisd()\" function needs to return positive values! (uisd(",theta,") <= 0)"))
      }
      return(result)
    }
    uisdFun <- Vectorize(uisdFun, "theta")
  }

  ####################
  # determine method:
  method <- match.arg(tolower(method), c("elir", "vr", "pr", "mtm.pt"))
  # switch:
  if (method == "elir") { #  "expected local-information-ratio (ELIR)" method:
    # specify function "i(p(theta))"
    # (equation (3) in Neuenschwander & al. (2020)):
    ip <- function(object, theta)
    {
      stopifnot(length(theta)==1, is.finite(theta))
      # compute Hessian:
      hessi <- hessian(function(x){object$dposterior(theta=x,
                                                     predict=TRUE,
                                                     log=TRUE)},
                       theta)
      # note the slight hack here:
      margin <- 6 * diff(object$summary[c("95% lower","95% upper"),"theta"])
      # (this should -roughly- correspond to a 20-sigma margin)
      if ((!is.finite(hessi))
          && (abs(theta-object$summary["median","theta"]) > margin)) {
        hessi <- 0.0
      }
      # (Hessian is set to zero in case of numerical problems
      #  AND a ~ 20-sigma difference from median)
      return(as.vector(-hessi))
    }
    ip <- Vectorize(ip, "theta")

    # compute expectation  (equation (7)):
    integrand <- function(x)
    {
      # NB: For numerical ease, the integrand is computed in a structured way.
      #     The integrand results as a product of 3 factors.
      #     Factors are computed one-by-one, _UNLESS_ one of the factors
      #     turned out as zero already.
      factors <- matrix(0.0, nrow=length(x), ncol=3,
                        dimnames=list(NULL, c("ip", "uisd2", "density")))
      # first, compute density:
      factors[,"density"] <- object$dposterior(theta=x, predict=TRUE)
      # second, compute i(p(theta)):
      zero <- (factors[,"density"] == 0.0)
      if (any(!zero)) {
        factors[!zero, "ip"] <- ip(object, theta=x[!zero])
      }
      # third, compute uisd^2:
      zero <- (factors[,"ip"] == 0.0)
      if (any(!zero)) {
        factors[!zero, "uisd2"] <- uisdFun(x[!zero])^2
      }
      # multiply 3 factors:
      result <- apply(factors, 1, prod)
      return(result)
    }
    expect <- integrate(integrand, lower=-Inf, upper=Inf)
    if (expect$message != "OK")
      warning(paste0("Problem computing expectation (\"", expect$message,"\")."))
    ESS <- expect$value
  } else if (method == "vr") {  # "variance ratio (VR)" method:
    # compute expectation  (equation (2)):
    if (is.vector(uisd)) {
      numerator <- uisd^2
    } else {
      integrand <- function(x)
      {
        # NB: For numerical ease, the integrand is computed in a structured way.
        #     The integrand results as a product of 2 factors (uisd and density).
        #     Factors are computed one-by-one, _UNLESS_ one of the factors
        #     turned out as zero already.
        factors <- matrix(0.0, nrow=length(x), ncol=2,
                          dimnames=list(NULL, c("uisd2", "density")))
        # first, compute density:
        factors[,"density"] <- object$dposterior(theta=x, predict=TRUE)
        # second, compute i(p(theta)):
        zero <- (factors[,"density"] == 0.0)
        if (any(!zero)) {
          factors[!zero, "uisd2"] <- uisdFun(x[!zero])^2
        }
        # multiply 2 factors:
        result <- apply(factors, 1, prod)
        return(result)
      }
      expect <- integrate(integrand, lower=-Inf, upper=Inf)
      if (expect$message != "OK")
        warning(paste0("Problem computing expectation (\"", expect$message,"\")."))
      numerator <- expect$value
    }
    denominator <- object$summary["sd","theta"]^2
    ESS <- numerator / denominator
  } else if (method == "pr") {  # "precision ratio (PR)" method:
    # compute expectation  (equation (2)):
    if (is.vector(uisd)) {
      denominator <- uisd^-2
    } else {
      integrand <- function(x)
      {
        # NB: For numerical ease, the integrand is computed in a structured way.
        #     The integrand results as a product of 2 factors (uisd and density).
        #     Factors are computed one-by-one, _UNLESS_ one of the factors
        #     turned out as zero already.
        factors <- matrix(0.0, nrow=length(x), ncol=2,
                          dimnames=list(NULL, c("uisd2", "density")))
        # first, compute density:
        factors[,"density"] <- object$dposterior(theta=x, predict=TRUE)
        # second, compute i(p(theta)):
        zero <- (factors[,"density"] == 0.0)
        if (any(!zero)) {
          factors[!zero, "uisd2"] <- uisdFun(x[!zero])^-2
        }
        # multiply 2 factors:
        result <- apply(factors, 1, prod)
        return(result)
      }
      expect <- integrate(integrand, lower=-Inf, upper=Inf)
      if (expect$message != "OK")
        warning(paste0("Problem computing expectation (\"", expect$message,"\")."))
      denominator <- expect$value
    }
    numerator <- object$summary["sd","theta"]^-2
    ESS <- numerator / denominator
  } else if (method == "mtm.pt"){  # "Morita-Thall-Mueller / Pennello-Thompson (MTM.PM)" method:
    denominator <- uisdFun(theta=object$summary["mode","theta"])^-2
    hessi <- hessian(function(x){object$dposterior(theta=x,
                                                   predict=TRUE,
                                                   log=TRUE)},
                     object$summary["mode","theta"])
    numerator <- as.numeric(-hessi)
    ESS <- numerator / denominator
  }
  return(ESS)
}


weightsplot <- function(x, ...)
{
  UseMethod("weightsplot")
}

weightsplot.bayesmeta <- function(x, individual=FALSE, ordered=TRUE,
                                  extramargin=4,
                                  priorlabel="prior mean", main, ...)
# Illustrate posterior mean weights (percentages)
# for overall mean or shrinkage estimates in a bar plot
# (See https://doi.org/10.1002/bimj.202000227 ).
# Numbers are taken from the bayesmeta object's
# "...$weights" or "...$weights.theta" elements.
{
  stopifnot(length(ordered)==1, is.logical(ordered),
            is.numeric(extramargin), length(extramargin)==1, all(is.finite(extramargin)),
            length(individual)==1)
  if (! (is.logical(individual) && (!individual))) { #  individual != FALSE
    indiv.logi <- TRUE       # (non-empty "individual" specification)
    if (is.numeric(individual))   indiv.which <- which(is.element(1:x$k, individual))
    if (is.character(individual)) indiv.which <- which(is.element(x$labels, match.arg(individual,x$labels)))
    if (length(indiv.which)==0) warning("cannot make sense of 'individual' argument: empty subset.")
  }
  else indiv.logi <- FALSE   # (the default (overall mean weight))
  mu.prior.proper <- all(is.finite(x$mu.prior))
  # create empty data frame to hold plotted data:
  plotdat <- cbind.data.frame("label"=rep(NA_character_, ifelse(mu.prior.proper, x$k+1, x$k)),
                              "weight"=NA_real_, "percentage"=NA_character_, stringsAsFactors=FALSE)
  if (!indiv.logi){ # fill in data for overall mean estimate
    plotdat[,"label"]  <- names(x$weights)
    plotdat[,"weight"] <- x$weights
    targetparam <- "overall mean estimate"
  } else {          # fill in data for shrinkage estimate
    plotdat[,"label"] <- rownames(x$weights.theta)
    plotdat[,"weight"] <- x$weights.theta[,indiv.which]
    targetparam <- paste0("shrinkage estimate \"", x$labels[indiv.which], "\"")
  }
  if (ordered) {    # sort rows in descending order
    plotdat[1:x$k,] <- plotdat[order(plotdat[1:x$k,"weight"], decreasing=TRUE),]
  }
  plotdat[,"percentage"] <- paste(sprintf("%.1f", plotdat[,"weight"] * 100), "%")
  maxweight <- max(plotdat[,"weight"]*100)
  if (missing(main)) {
    main <- paste0("posterior mean weights (",targetparam,")")
  }
  # set figure margins:
  if (extramargin != 0) {
    parmar <- par("mar")
    on.exit(par(mar=parmar))
    par(mar = parmar + c(0, extramargin, 0, 0))
  }
  bp <- graphics::barplot(rev(plotdat[,"weight"])*100, horiz=TRUE,
                          names.arg=rev(plotdat[,"label"]), las=1,
                          xlim=c(0, maxweight * 1.15),
                          xlab="weight (%)", ylab="", main=main, ...)
  graphics::text(rev(plotdat[,"weight"])*100 + maxweight*0.02, bp[,1],
                 rev(plotdat[,"percentage"]), adj=c(0,0.5))
  invisible(plotdat)
}


traceplot <- function(x, ...)
{
  UseMethod("traceplot")
}


traceplot.bayesmeta <- function(x, mulim, taulim, ci=FALSE,
                                rightmargin=8, col=rainbow(x$k), ...)
{
  stopifnot(missing(mulim) || (length(mulim) == 2),
            missing(taulim) || (length(taulim) <= 2),
            rightmargin >= 0, length(col) == x$k)
  # convert "taulim" and "mulim" arguments
  # to eventual "taurange" and "murange" vectors:
  if (!missing(taulim) && all(is.finite(taulim))) {
    if ((length(taulim)==2) && (taulim[1]>=0) && (taulim[2]>taulim[1]))
      taurange <- taulim
    else if ((length(taulim)==1) && (taulim>0))
      taurange <- c(0, taulim)
    else
      taurange <- c(0, x$qposterior(tau=0.995)*1.1)
  } else {
    taurange <- c(0, x$qposterior(tau=0.995)*1.1)
  }
  
  if (!missing(mulim) && (all(is.finite(mulim)) && (mulim[1] < mulim[2]))) {
    murange <- mulim
  } else {
    cm <- x$cond.moment(tau=taurange[2], indiv=TRUE)
    if (ci){
      murange <- range(c(cm[,"mean",]-q975*cm[,"sd",], cm[,"mean",]+q975*cm[,"sd",]))
    } else {
      murange <- range(cm[,"mean",])
    }
    murange <- murange + c(-1,1)*diff(murange)*0.05
  }
  
  vertlines <- pretty(taurange)
  gridcol <- "grey85"
  q975 <- qnorm(0.975)

  mutrace <- function(x)
  {
    # range of tau values:
    tau <- seq(max(c(0,taurange[1]-0.1*diff(taurange))),
                   taurange[2]+0.1*diff(taurange), le=200)
    cm.overall <- x$cond.moment(tau=tau)
    cm.indiv   <- x$cond.moment(tau=tau, indiv=TRUE)
    plot(taurange, murange,         
         type="n", axes=FALSE, xlab="", ylab="effect", main="", ...)
    abline(v=vertlines, col=gridcol)
    abline(h=pretty(murange), col=gridcol)
    abline(v=0, col=grey(0.40))
    # grey shading:
    if (ci) {
      for (i in 1:x$k) {
        polygon(c(tau, rev(tau)),
                c(cm.indiv[,"mean",i] - q975*cm.indiv[,"sd",i],
                  rev(cm.indiv[,"mean",i] + q975*cm.indiv[,"sd",i])),
                col=grey(0.75, alpha=0.25), border=NA)
      }
      polygon(c(tau, rev(tau)),
              c(cm.overall[,"mean"] - q975*cm.overall[,"sd"],
                rev(cm.overall[,"mean"] + q975*cm.overall[,"sd"])),
              col=grey(0.75, alpha=0.25), border=NA)
    }    
    # individual estimates:
    matlines(tau, cm.indiv[,"mean",], col=col, lty=1)
    if (ci) {
      matlines(tau, cm.indiv[,"mean",]-q975*cm.indiv[,"sd",], col=col, lty=3)
      matlines(tau, cm.indiv[,"mean",]+q975*cm.indiv[,"sd",], col=col, lty=3)
    }
    # overall mean:
    lines(tau, cm.overall[,"mean"], col="black", lty=2, lwd=1.5)
    if (ci) {
      lines(tau, cm.overall[,"mean"]-q975*cm.overall[,"sd"], col="black", lty=3, lwd=1.5)
      lines(tau, cm.overall[,"mean"]+q975*cm.overall[,"sd"], col="black", lty=3, lwd=1.5)
    }
    axis(2)
    for (i in 1:x$k)
      axis(side=4, at=cm.indiv[length(tau),"mean",i],
           labels=x$labels[i], tick=FALSE,
           col.axis=col[i], las=1)
    axis(side=4, at=cm.overall[length(tau),"mean"],
         labels="overall mean", tick=FALSE, las=1)
    invisible()
  }
  
  taumarginal <- function(x)
  # NB: function is (essentially) identical to the one within "plot.bayesmeta()"
  {
    # range of tau values:
    tau <- seq(max(c(0,taurange[1]-0.1*diff(taurange))),
                   taurange[2]+0.1*diff(taurange), le=200)
    # corresponding posterior density:
    dens <- x$dposterior(tau=tau)
    # empty plot:
    maxdens <- max(dens[is.finite(dens)],na.rm=TRUE)
    plot(c(taurange[1],taurange[2]), c(0,maxdens),         
         type="n", axes=FALSE, xlab="", ylab="", main="")
    abline(v=vertlines, col=gridcol)
    # "fix" diverging density:
    dens[!is.finite(dens)] <- 10*maxdens
    # light grey shaded contour for density across whole range:
    polygon(c(0,tau,max(tau)), c(0,dens,0), border=NA, col=grey(0.90))
    # dark grey shaded contour for density within 95% bounds:
    indi <- ((tau>=x$summary["95% lower","tau"]) & (tau<=x$summary["95% upper","tau"]))
    polygon(c(rep(x$summary["95% lower","tau"],2), tau[indi], rep(x$summary["95% upper","tau"],2)),
            c(0, min(c(x$dposterior(tau=x$summary["95% lower","tau"]), 10*maxdens)),
              dens[indi], x$dposterior(tau=x$summary["95% upper","tau"]), 0),
            border=NA, col=grey(0.80))
    # vertical line at posterior median:
    lines(rep(x$summary["median","tau"],2), c(0,x$dposterior(tau=x$summary["median","tau"])), col=grey(0.6))
    # actual density line:
    lines(tau, dens, col="black")
    # x-axis, y-axis:
    abline(h=0, v=0, col=grey(0.40))
    # add axes, labels, bounding box, ...
    mtext(side=1, line=par("mgp")[1], expression("heterogeneity "*tau))
    #mtext(side=2, line=par("mgp")[2], expression("marginal posterior density"))
    axis(1)#; box()
    invisible()
  }

  # make sure to re-set graphical parameters later:
  prevpar <- par(no.readonly=TRUE)
  on.exit(par(prevpar))
  # generate actual plot:
  graphics::layout(rbind(1,2), heights=c(2,1))
  par(mar=c(-0.1,3,0,rightmargin)+0.1, mgp=c(2.0, 0.8, 0))
  mutrace(x)
  par(mar=c(3,3,-0.1,rightmargin)+0.1)
  taumarginal(x)
  graphics::layout(1)
  par(mar=c(5,4,4,2)+0.1)
  invisible()
}
