#
#    bayesmeta, an R package for Bayesian random-effects meta-analysis.
#    Copyright (C) 2023  Christian Roever
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


kldiv <- function(mu1, mu2, sigma1, sigma2, symmetrized=FALSE)
# Kullback-Leibler divergence (relative entropy) of two normal distributions
#   I(1:2)
# or /symmetrized/ divergence
#   J(1,2) = I(1:2) + I(2:1)
#
# Arguments:
#   mu1, mu2       :  the two distributions' mean vectors 
#   sigma1, sigma2 :  the two distributions' covariance matrices
#   symmetrized    :  indicator of whether to compute symmetrized divergence
#
{
  stopifnot(is.vector(mu1), is.vector(mu2),
            length(mu1) == length(mu2),
            is.numeric(mu1), is.numeric(mu2),
            all(is.finite(mu1)), all(is.finite(mu2)),
            is.matrix(sigma1), is.matrix(sigma2),
            nrow(sigma1) == ncol(sigma1),
            all(dim(sigma1) == dim(sigma2)))
  k <- length(mu1)
  sigma2inv <- solve(sigma2)
  KL <- 0.5 * (sum(diag(sigma2inv %*% sigma1))
               + (t(mu1-mu2) %*% sigma2inv %*% (mu1-mu2)) - k)
  if (symmetrized) {
    sigma1inv <- solve(sigma1)
    KL <- KL + 0.5 * (sum(diag(sigma1inv %*% sigma2))
                      + (t(mu2-mu1) %*% sigma1inv %*% (mu2-mu1)) - k)
  } else {
    #KL <- KL + 0.5 * (log(det(sigma2)) - log(det(sigma1)))
    KL <- KL + 0.5 * (determinant(sigma2, logarithm=TRUE)$modulus
                      - determinant(sigma1, logarithm=TRUE)$modulus)
    # (NB: this term cancels out for the symmetrized variant)
  }
  return(as.vector(KL))
}


################################################################################

beta.convert <- function(beta, which.beta, d, betanames)
# checks "beta" and "which.beta" arguments for possible errors,
# and returns a d-column "beta" _matrix_
# based on a range of possible input specifications.
{
  stopifnot(is.matrix(beta) | is.data.frame(beta) | is.vector(beta))
  if (missing(which.beta)) { # default "which.beta"
    betaIdx <- 1:d
  } else {  # try to make sense of "which.beta" argument; general sanity checks:
    #if (!is.vector(which.beta)
    #    || (!is.element(class(which.beta),
    #                    c("numeric", "integer", "logical", "character"))
    #        | (length(which.beta) > d))) {
    #  warning("Cannot make sense of \"which.beta\" argument (1).")
    #}
    if (!is.vector(which.beta)
        || (!inherits(which.beta,
                      c("numeric", "integer", "logical", "character"))
            | (length(which.beta) > d))) {
      warning("Cannot make sense of \"which.beta\" argument (1).")
    }
     if (is.character(which.beta)){       # character "which.beta"
      if (!all(is.element(which.beta, betanames))
          | any(duplicated(which.beta)))
        warning("Cannot make sense of \"which.beta\" argument (2).")
      betaIdx <- rep(NA_integer_, length(which.beta))
      for (i in 1:length(which.beta))
        betaIdx[i] <- which(is.element(betanames, which.beta[i]))
    } else if (is.logical(which.beta)) { # logical "which.beta"
      if ((length(which.beta) != d) | any(!is.finite(which.beta))) 
        warning("Cannot make sense of \"which.beta\" argument (3).")
      betaIdx <- which(which.beta)
    } else {                             # numeric/integer "which.beta"
      if (!all(is.element(which.beta, 1:d))
          | any(duplicated(which.beta)))
        warning("Cannot make sense of \"which.beta\" argument (4).")
      betaIdx <- which.beta
    }
  }
  # now have a numeric "betaIdx" index vector, indicating which exact variable
  # the columns (or elements) of the supplied "beta" argument correspond to.
  if (is.data.frame(beta)) beta <- as.matrix(beta)
  if (is.vector(beta)) {
    if (length(betaIdx)==1)
      beta <- matrix(beta, ncol=1)
    else if (length(beta) == length(betaIdx))
      beta <- matrix(beta, nrow=1)
    else
      warning("non-conformable \"beta\" argument.")
  }
  stopifnot(is.numeric(beta))
  # ensure either all-NA or all-finite values per column:
  stopifnot(all(apply(beta, 2, function(x){all(is.na(x)) | all(is.finite(x))})))
  # generate "new" beta matrix of standardized format:
  betaMatrix <- matrix(NA_real_, ncol=d, nrow=nrow(beta),
                       dimnames=list(NULL, betanames))
  for (i in 1:length(betaIdx))
    betaMatrix[,betaIdx[i]] <- beta[,i]
  return(betaMatrix)
}


################################################################################


bmr <- function(y,...)
{
  UseMethod("bmr")
}


bmr.default <- function(y, sigma, labels = names(y),
                        X = matrix(1.0, nrow=length(y), ncol=1,
                                   dimnames=list(labels,"intercept")),
                        tau.prior = "uniform",
                        beta.prior.mean = NULL,
                        beta.prior.sd   = NULL,
                        beta.prior.cov  = diag(beta.prior.sd^2,
                                               nrow=length(beta.prior.sd),
                                               ncol=length(beta.prior.sd)),
                        interval.type = c("shortest", "central"),
                        delta = 0.01, epsilon = 0.0001,
                        rel.tol.integrate=2^16*.Machine$double.eps,
                        abs.tol.integrate=0.0,
                        tol.uniroot=rel.tol.integrate, ...)
{
  ptm <- proc.time()
  if (is.vector(X)) {
    X <- matrix(X, ncol=1)
  }
  if (!is.vector(y))     y     <- as.vector(y)
  if (!is.vector(sigma)) sigma <- as.vector(sigma)
  stopifnot(is.vector(y), is.numeric(y), all(is.finite(y)),
            is.vector(sigma), is.numeric(sigma), all(is.finite(sigma)), all(sigma>0),
            is.matrix(X), is.numeric(X), all(is.finite(X)),
            length(sigma)==length(y), nrow(X)==length(y),
            (is.function(tau.prior) | (is.character(tau.prior) && (length(tau.prior)==1))))
  if (!requireNamespace("mvtnorm", quietly=TRUE))
    stop("required 'mvtnorm' package not available!")
  if (ncol(X) > nrow(X)) {
    warning(paste0("input matrix \"X\" has more columns than rows ",
                   "(ncol(X)=", ncol(X), " > ", nrow(X), "=nrow(X))."))
  }
  Xrank <- qr(X)$rank 
  if (Xrank < ncol(X)) {
    warning(paste0("input matrix \"X\" is not of full rank ",
                   "(rank(X)=",Xrank," < ",ncol(X),"=ncol(X))."))
  }
  if (is.null(colnames(X))) { # generate column names, if not provided:
    colnames(X) <- sprintf("var%03d", 1:ncol(X))
  } else {                    # otherwise ensure proper naming:
    colnames(X) <- make.names(colnames(X), unique=TRUE)
  }
  stopifnot(!is.element("tau", colnames(X)))
  betanames <- colnames(X)
  k <- length(y) # number of studies
  d <- ncol(X)   # number of covariables
  if (is.null(labels))
    labels <- sprintf("%02d", 1:k)
  if (is.null(names(y)))
    names(y) <- labels
  if (is.null(names(sigma)))
    names(sigma) <- labels
  if (is.null(rownames(X)))
    rownames(X) <- labels
  # specify extended "y", "sigma" and "X" vectors/matrices
  # INCLUDING prior information:
  if (!is.null(beta.prior.mean)) {
    # sanity checks:
    stopifnot(length(beta.prior.mean) == d)
    if (!is.null(beta.prior.sd)) {
      stopifnot(length(beta.prior.sd) == d,
                all(beta.prior.sd > 0))
      names(beta.prior.sd) <- betanames
    }
    stopifnot(length(beta.prior.mean) == ncol(X),
              is.matrix(beta.prior.cov),
              nrow(beta.prior.cov) == length(beta.prior.mean),
              ncol(beta.prior.cov) == nrow(beta.prior.cov))
    names(beta.prior.mean) <- betanames
    rownames(beta.prior.cov) <- colnames(beta.prior.cov) <- betanames
    # determine finite prior variances
    # ( == subset of beta parameters with prior information) 
    beta.prior.proper <- is.finite(diag(beta.prior.cov))
    if (!any(beta.prior.proper))
      warning("need to specify at least one finite beta prior variance parameter!")
    if (any(!is.finite(beta.prior.mean[beta.prior.proper])))
      warning("need to specify finite beta prior mean parameter values!")
    yp     <- beta.prior.mean[beta.prior.proper]
    Xp     <- diag(d)[beta.prior.proper, beta.prior.proper]
    Sigmap <- beta.prior.cov[beta.prior.proper,beta.prior.proper]
  } else {  # (improper uniform prior)
    yp <- Xp <- Sigmap <- NULL
    beta.prior.proper <- rep(FALSE, d)
  }

  # assemble prior parameters (also for output later)
  names(beta.prior.proper) <- betanames
  betaprior <- list("mean"       =rep(NA_real_, d),
                    "covariance" =diag(Inf, nrow=d, ncol=d))
  names(betaprior$mean) <- betanames
  rownames(betaprior$covariance) <- colnames(betaprior$covariance) <- betanames
  if (any(beta.prior.proper)) {
    betaprior$mean[beta.prior.proper] <- beta.prior.mean[beta.prior.proper]
    betaprior$covariance[beta.prior.proper,beta.prior.proper] <- beta.prior.cov[beta.prior.proper,beta.prior.proper]
  }
  
  interval.type <- match.arg(interval.type)

  sigma2hat <- (k-1)*sum(1/sigma^2) / (sum(1/sigma^2)^2 - sum(1/sigma^4)) # Higgins/Thompson (2002), eqn. (9)
  
  # check for extreme sigma values:
  maxratio <- max(sigma) / min(sigma)
  if (maxratio > 1000)
    warning(paste0("Ratio of largest over smallest standard error (sigma) is ", sprintf("%.0f",maxratio), ". Extreme values may lead to computational problems."))

  # digest "tau.prior" argument:
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
  
  dprior <- function(tau, beta, which.beta, log=FALSE)
  {
    stopifnot(!missing(tau) | !missing(beta))
    # beta density (if required)
    if (!missing(beta)) {  # marginal beta density
      betaMatrix <- beta.convert(beta, which.beta, d, betanames)
      if (!missing(tau))
        stopifnot(length(tau) == nrow(betaMatrix))
      beta.provided <- apply(betaMatrix, 2, function(b){all(is.finite(b))})
      relevant <- (beta.prior.proper & beta.provided)
      if (any(relevant))
        logbetadens <- mvtnorm::dmvnorm(x=betaMatrix[, relevant, drop=FALSE],
                                        mean=betaprior$mean[relevant],
                                        sigma=betaprior$covariance[relevant,relevant,drop=FALSE],
                                        log=TRUE)
      else
        logbetadens <- rep(0.0, nrow(betaMatrix))
    } else {                    # joint density
      logbetadens <- rep(0.0, length(tau))
    }
    # tau density (if required)
    if (!missing(tau)) {   # marginal tau density
      if (is.element("log", names(formals(tau.prior))))
        logtaudens <- apply(matrix(tau,ncol=1), 1, tau.prior, log=TRUE) - log(tau.prior.integral)
      else
        logtaudens <- log(apply(matrix(tau,ncol=1), 1, tau.prior)) - log(tau.prior.integral)
    } else {
      logtaudens <- rep(0.0, nrow(betaMatrix))
    }
    # combine:
    logpriordens <- logtaudens + logbetadens
    if (log) result <- logpriordens
    else result <- exp(logpriordens)
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

  tau.posterior.integral  <- 1.0 

  dposterior <- function(tau, beta, which.beta, log=FALSE)
  # posterior density for tau and/or beta
  {
    stopifnot(missing(tau)  || (is.vector(tau) & is.numeric(tau)),
              missing(beta) || (is.matrix(beta) | is.data.frame(beta) | is.vector(beta)))
    if (!missing(beta)) {
      betaMatrix <- beta.convert(beta, which.beta, d, betanames)
      stopifnot(missing(tau) || (length(tau) == nrow(betaMatrix)))
    }
    if ((!missing(tau) & !missing(beta))
        && (nrow(betaMatrix) != length(tau)))
      warning("non-conformable arguments \"tau\" and \"beta\".")
    
    taulogpostdens <- function(t)
    # marginal posterior density of tau  (marginalized over beta)
    {
      stopifnot(length(t) == 1)
      logdens <- -Inf
      if (is.finite(t) && ((t >= 0) & (t < Inf))) {
        if (!is.null(Sigmap)) { # set up block diagonal matrices:
          sigmaTau    <- cbind(rbind(diag(sigma^2 + t^2, nrow=k, ncol=k),
                                     matrix(0.0, nrow=length(yp), ncol=k)),
                               rbind(matrix(0.0, nrow=k, ncol=length(yp)),
                                     Sigmap))
          sigmaTauInv <- cbind(rbind(diag(1 / (sigma^2 + t^2), nrow=k, ncol=k),
                                     matrix(0.0, nrow=length(yp), ncol=k)),
                               rbind(matrix(0.0, nrow=k, ncol=length(yp)),
                                     solve(Sigmap)))
        } else {
          sigmaTau    <- diag(sigma^2 + t^2, nrow=k, ncol=k)
          sigmaTauInv <- diag(1 / (sigma^2 + t^2), nrow=k, ncol=k)
        }
        Vbeta       <- solve(t(rbind(X,Xp)) %*% sigmaTauInv %*% rbind(X,Xp))
        betaHat     <- Vbeta %*% t(rbind(X,Xp)) %*% sigmaTauInv %*% c(y,yp)
        residual    <- c(y,yp) - rbind(X,Xp) %*% betaHat
        logdens  <- (0.5 * ((d-k) * log(2*pi)
                            #-log(det(sigmaTau)) + log(det(Vbeta))
                            -determinant(sigmaTau, logarithm=TRUE)$modulus
                            +determinant(Vbeta, logarithm=TRUE)$modulus
                            -(t(residual) %*% sigmaTauInv %*% residual))
                     + log(tau.prior(t)))
      }
      return(logdens - log(tau.posterior.integral))
    }

    betalogpostdens <- function(tau, beta)
    # joint posterior density of tau and (some) beta parameters
    # (potentially also marginalized over tau via grid summation)
    {
      idx <- (!is.na(beta)) # only consider those variables with non-NA entries.
      if (missing(tau)) { # marginal density of beta (marginalized over tau)
        logdens <- numeric(length(support$weight))
        for (i in 1:length(support$weight)) {
          # NEED TO SPEED UP HERE
          # (need to dissect tau conditional and tau marginal into separate functions).
          covmat <- matrix(support$covariance[i,,], nrow=d, ncol=d)
          logdens[i] <- mvtnorm::dmvnorm(x=beta[idx], mean=support$mean[i,idx],
                                         sigma=covmat[idx,idx,drop=FALSE], log=TRUE)
        }
        logdens <- log(sum(exp(log(support$weight) + logdens)))
      } else {            # joint density of (tau,beta)
        if (tau>=0) {
          cm <- conditionalmoments(tau=tau)
          logdens <- (taulogpostdens(tau)
                      + mvtnorm::dmvnorm(x=beta[idx], mean=cm$mean[idx],
                                         sigma=cm$covariance[idx,idx,drop=FALSE], log=TRUE))
        } else {
          logdens <- -Inf
        }
      }
      return(logdens)
    }

    if (missing(beta)){        # marginal posterior of tau:
      ldens <- apply(matrix(tau, ncol=1), 1, taulogpostdens)
    } else if (missing(tau)) { # marginal posterior of beta
      ldens <- apply(betaMatrix, 1, function(x){betalogpostdens(beta=x)})
    } else {                   # joint posterior of (tau, beta)
      ldens <- apply(cbind(tau, betaMatrix), 1,
                     function(x){betalogpostdens(tau=x[1], beta=x[-1])})
    }
    
    if (log) result <- ldens
    else     result <- exp(ldens)
    return(result)
  }

  # determine tau posterior's normalizing constant:
  tpint <- integrate(function(t){return(dposterior(tau=t))},
                     lower=0, upper=Inf,
                     rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate,
                     stop.on.error=FALSE)
  if (tpint$message != "OK")
    warning("Failed integrating tau posterior density (\"", tpint$message, "\")!")
  tau.posterior.integral <- tpint$value

  pposterior <- function(tau, beta, which.beta)
  # posterior cumulative distribution function (CDF)
  # for tau or beta
  {
    stopifnot(xor(missing(tau), missing(beta)),
              missing(tau)  || (is.vector(tau) & is.numeric(tau)),
              missing(beta) || (is.matrix(beta) | is.data.frame(beta) | is.vector(beta)))
    if (!missing(beta)) {
      if (!missing(which.beta)) { # general sanity checks
        #if (!is.vector(which.beta)
        #    || (!is.element(class(which.beta),
        #                    c("numeric", "integer", "logical", "character"))
        #        | ((class(which.beta)=="logical")
        #            && ((length(which.beta) != d) | (sum(which.beta)!=1)))
        #        | ((is.element(class(which.beta),
        #                       c("numeric", "integer", "character"))
        #            & (length(which.beta) != 1))))) {
        #  warning("Cannot make sense of \"which.beta\" argument (1).")
        #}
        if (!is.vector(which.beta)
            || (!inherits(which.beta,
                          c("numeric", "integer", "logical", "character"))
                | (inherits(which.beta, "logical")
                    && ((length(which.beta) != d) | (sum(which.beta)!=1)))
                | ((inherits(which.beta,
                             c("numeric", "integer", "character"))
                    & (length(which.beta) != 1))))) {
          warning("Cannot make sense of \"which.beta\" argument (1).")
        }
        if (is.character(which.beta)){       # character "which.beta"
          if (!is.element(which.beta, betanames))
            warning("Cannot make sense of \"which.beta\" argument (2).")
          betaIdx <- which(is.element(betanames, which.beta))
        } else if (is.logical(which.beta)) { # logical "which.beta"
          if (any(!is.finite(which.beta)))
            warning("Cannot make sense of \"which.beta\" argument (3).")
          betaIdx <- which(which.beta)
        } else {                             # numeric/integer "which.beta"
          if (!(is.element(which.beta, 1:d)))
            warning("Cannot make sense of \"which.beta\" argument (4).")
          betaIdx <- which.beta
        }
      }
      # now have a numeric "betaIdx" index vector, indicating which exact variable
      # the supplied "beta" argument corresponds to.
      if (is.data.frame(beta)) beta <- as.matrix(beta)
      if (is.vector(beta)) {
        beta <- matrix(beta, ncol=1)
      }
      stopifnot(is.numeric(beta))
      # generate "new" beta matrix of standardized format:
      if (!missing(which.beta)) {
        betaMatrix <- matrix(NA_real_, ncol=d, nrow=nrow(beta),
                             dimnames=list(NULL, betanames))
        betaMatrix[,betaIdx] <- beta[,1]
      } else {
        stopifnot(ncol(beta) == d)
        betaMatrix <- beta
      }
    }
    stopifnot(xor(missing(tau), missing(beta)),
              missing(tau)  || (is.vector(tau) & is.numeric(tau)))
    
    taucdf <- function(t)
    {
      prob <- NA_real_
      if (t <= 0) {
        prob <- 0.0
      } else if (t==Inf) {
        prob <- 1.0
      } else if (is.finite(t)) {
        prob <- integrate(function(x){dposterior(tau=x)}, lower=0, upper=t,
                          rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value
      }
      return(prob)
    }
    
    betacdf <- function(b)
    {
      idx <- (!is.na(b))
      prob <- sum(exp(pnorm(b[idx], mean=support$mean[,idx],
                            sd=sqrt(support$covariance[,idx,idx]), log.p=TRUE)
                      + log(support$weight)))
      return(prob)
    }
    
    if (!missing(tau)) {
      cdf <- apply(matrix(tau, ncol=1), 1, taucdf)
    } else {
      finitebeta <- rowSums(!is.na(betaMatrix))
      if (any(finitebeta != 1)) {
        warning("need to specify exactly 1 finite beta value per row")
      }
      cdf <- apply(betaMatrix, 1, betacdf)
    }
    return(cdf)
  }

  qposterior <- function(tau.p, beta.p, which.beta)
  # posterior quantile function (inverse CDF) for tau or beta
  {
    stopifnot(xor(missing(tau.p), missing(beta.p)),
              missing(tau.p)  || (is.vector(tau.p) & is.numeric(tau.p)),
              missing(beta.p) || (is.matrix(beta.p) | is.data.frame(beta.p) | is.vector(beta.p)))
    if (!missing(beta.p)) {
      if (!missing(which.beta)) { # general sanity checks
        #if (!is.vector(which.beta)
        #    || (!is.element(class(which.beta),
        #                    c("numeric", "integer", "logical", "character"))
        #        | ((class(which.beta)=="logical")
        #            & ((length(which.beta) != d) | (sum(which.beta) != 1)))
        #        | ((is.element(class(which.beta),
        #                       c("numeric", "integer", "character"))
        #            & (length(which.beta) != 1))))) {
        #  warning("Cannot make sense of \"which.beta\" argument (1).")
        #}
        if (!is.vector(which.beta)
            || (!inherits(which.beta,
                          c("numeric", "integer", "logical", "character"))
                | (inherits(which.beta,"logical")
                   & ((length(which.beta) != d) | (sum(which.beta) != 1)))
                | ((inherits(which.beta,
                             c("numeric", "integer", "character"))
                    & (length(which.beta) != 1))))) {
          warning("Cannot make sense of \"which.beta\" argument (1).")
        }
        if (is.character(which.beta)){       # character "which.beta"
          if (!is.element(which.beta, betanames))
            warning("Cannot make sense of \"which.beta\" argument (2).")
          betaIdx <- which(is.element(betanames, which.beta))
        } else if (is.logical(which.beta)) { # logical "which.beta"
          if (any(!is.finite(which.beta)))
            warning("Cannot make sense of \"which.beta\" argument (3).")
          betaIdx <- which(which.beta)
        } else {                             # numeric/integer "which.beta"
          if (!(is.element(which.beta, 1:d)))
            warning("Cannot make sense of \"which.beta\" argument (4).")
          betaIdx <- which.beta
        }
      }
      # now have a numeric "betaIdx" index vector, indicating which exact variable
      # the supplied "beta.p" argument corresponds to.
      if (is.data.frame(beta.p)) beta.p <- as.matrix(beta.p)
      if (is.vector(beta.p)) {
        beta.p <- matrix(beta.p, ncol=1)
      }
      stopifnot(is.numeric(beta.p))
      # generate "new" beta matrix of standardized format:
      if (!missing(which.beta)) {
        stopifnot(ncol(beta.p) == 1)
        betaMatrix <- matrix(NA_real_, ncol=d, nrow=nrow(beta.p),
                             dimnames=list(NULL, betanames))
        betaMatrix[,betaIdx] <- beta.p[,1]
      } else {
        stopifnot(ncol(beta.p) == d)
        betaMatrix <- beta.p
      }
    }
    
    tauquant <- function(p)
    {
      tau <- NA_real_
      if (is.finite(p)) {
        if (p==0) {
          tau <- 0.0
        } else if (p==1) {
          tau <- Inf
        } else if ((p>0) & (p<1)) {
          upper <- 1.0
          while (pposterior(tau=upper) < p) upper <- upper * 2
          tau <- uniroot(function(xx){return(pposterior(tau=xx)-p)},
                         lower=0, upper=upper, tol=tol.uniroot)$root
        }
      }
      return(tau)
    }
    
    betaquant <- function(p, wb)
    {
      idx <- (!is.na(p))
      beta <- NA_real_
      if (is.finite(p[idx])) {
        if (p[idx]==0) {
          beta <- -Inf
        } else if (p[idx]==1) {
          beta <- Inf
        } else if ((p[idx]>0) & (p[idx]<1)) {
          #lower <- upper <- rep(NA,d)
          #lower[idx] <- upper[idx] <- median(support$mean[,idx])
          #step <- sqrt(median(support$covariance[,idx,idx]))
          #lower[idx] <- lower[idx] - 2*step
          #upper[idx] <- upper[idx] + 2*step
          #while (pposterior(beta=lower) > p[idx]) lower[idx] <- lower[idx] - step
          #while (pposterior(beta=upper) < p[idx]) upper[idx] <- upper[idx] + step
          #auxfun <- function(xx) {
          #  bvec      <- rep(NA_real_, d)
          #  bvec[idx] <- xx
          #  return(pposterior(beta=bvec)-p[idx])
          #}
          #beta <- uniroot(auxfun, lower=lower[idx], upper=upper[idx],
          #                tol=tol.uniroot)$root
          lower <- upper <- stats::median(support$mean[,idx])
          step <- sqrt(stats::median(support$covariance[,idx,idx]))
          lower <- lower - 2*step
          upper <- upper + 2*step
          while (pposterior(beta=lower, which=idx) > p[idx]) lower <- lower - step
          while (pposterior(beta=upper, which=idx) < p[idx]) upper <- upper + step
          auxfun <- function(xx) {
            bvec      <- rep(NA_real_, d)
            bvec[idx] <- xx
            return(pposterior(beta=bvec)-p)
          }
          beta <- uniroot(function(xx){pposterior(beta=xx, which=idx)-p[idx]},
                          lower=lower, upper=upper, tol=tol.uniroot)$root
        }
      }
      return(beta)
    }

    if (!missing(tau.p)) {
      result <- apply(matrix(tau.p, ncol=1), 1, tauquant)
    } else {
      finitebetap <- rowSums(!is.na(beta.p))
      if (any(finitebetap != 1)) {
        warning("need to specify exactly one finite beta value per row")
      }
      result <- apply(betaMatrix, 1, betaquant)      
    }
    return(result)
  }

  post.interval <- function(tau.level, beta.level, which.beta,
                            method=interval.type)
  {
    stopifnot(xor(missing(tau.level), missing(beta.level)),
              missing(tau.level)  || (is.numeric(tau.level) & (length(tau.level)==1)
                                      & (tau.level>0) & (tau.level<1)),
              missing(beta.level) || (is.numeric(beta.level) & (length(beta.level)==1)
                                      & (beta.level>0) & (beta.level<1)
                                      & !missing(which.beta)))
    if (!missing(tau.level)) { # credible interval for tau
      if (tau.level < 0.5) warning("\"tau.level\" below 50% ... intended?")
      if (method=="central") { # "central", equal-tailed interval:
        result <- qposterior(tau.p = c((1-tau.level)/2, 1-(1-tau.level)/2))
      } else {                 # "shortest" interval:
        intwidth <- function(left)
        {
          pleft <- pposterior(tau = left)
          right <- qposterior(tau.p = tau.level + pleft)
          return(right-left)
        }
        opti <- optimize(intwidth, lower=0, upper=qposterior(tau.p=1-tau.level))$minimum
        # catch marginal minimum:
        if (intwidth(0) < intwidth(opti))
          result <- c(0, qposterior(tau=tau.level))
        else
          result <- c(opti, qposterior(tau=tau.level+pposterior(tau=opti)))
      }
    } else {                   # credible interval for beta
      if (beta.level < 0.5) warning("\"beta.level\" below 50% ... intended?")
      if (method=="central") { # "central", equal-tailed interval:
        result <- qposterior(beta.p = c((1-beta.level)/2, 1-(1-beta.level)/2),
                             which.beta = which.beta)
      } else {                 # "shortest" interval:
        intwidth <- function(left)
        {
          pleft <- pposterior(beta = left, which.beta = which.beta)
          right <- qposterior(beta.p = beta.level + pleft, which.beta = which.beta)
          return(right-left)
        }
        opti <- optimize(intwidth,
                         lower=qposterior(beta.p=(1-beta.level)/50,
                                          which.beta = which.beta),
                         upper=qposterior(beta.p=1-beta.level,
                                          which.beta = which.beta))$minimum
        result <- c(opti, qposterior(beta.p = beta.level+pposterior(beta=opti,which.beta = which.beta),
                                     which.beta = which.beta))
      }
    }    
    attr(result, "interval.type") <- method
    return(result)
  }
  
  rposterior <- function(n=1, tau.sample=TRUE)
  {
    stopifnot(n>=1, n==round(n))
    if (tau.sample) {
      result <- matrix(NA_real_, nrow=n, ncol=d+1,
                       dimnames=list(NULL, c("tau", betanames)))
      # draw from tau's marginal posterior:
      result[,"tau"] <- qposterior(tau.p = runif(n=n, min=0, max=1))
      # draw beta values conditional on tau:
      rcondbeta <- function(tau)
      {
        cm <- conditionalmoments(tau)
        return(mvtnorm::rmvnorm(1, mean=cm$mean, sigma=cm$covariance))
      }
      result[, betanames] <- t(apply(result[,"tau",drop=FALSE], 1, rcondbeta))
    } else {
      # (hopefully) faster sampling schedule (based on DIRECT support grid)
      # in case tau samples are not necessary:
      result <- matrix(NA_real_, nrow=0, ncol=d,
                       dimnames=list(NULL, betanames))
      # determine numbers of samples from each mixture component:
      freq <- as.vector(stats::rmultinom(1, size=n, prob=support$weight))
      # draw actual (mv-normal) samples component-wise:
      for (i in 1:length(support$weight)) {
        if (freq[i] > 0) {
          covmat <- matrix(support$covariance[i,,], nrow=d, ncol=d)
          result <- rbind(result,
                          mvtnorm::rmvnorm(freq[i], mean=support$mean[i,],
                                          sigma=covmat))
        }
      }
      # re-shuffle rows:
      result <- result[sample(x=n, size=n, replace=FALSE),]
    }
    return(result)
  }

  dpredict <- function(theta, x, mean=TRUE, log=FALSE)
  {
    stopifnot(!missing(theta), is.vector(theta), is.numeric(theta), all(is.finite(theta)),
              !missing(x), is.vector(x), is.numeric(x),
              length(x) == d, all(is.finite(x)))
    predMom <- pred.moments(x=x, mean=mean)
    result <- rep(NA_real_, length(theta))
    logwts <- log(support$weight)
    for (i in 1:length(theta))
      result[i] <- sum(exp(logwts + dnorm(x=theta[i], mean=predMom[,"mean"], sd=predMom[,"sd"], log=TRUE)))
    if (log) result <- log(result)
    return(result)
  }

  ppredict <- function(theta, x, mean=TRUE)
  {
    stopifnot(!missing(theta), is.vector(theta), is.numeric(theta), all(is.finite(theta)),
              !missing(x), is.vector(x), is.numeric(x),
              length(x) == d, all(is.finite(x)))
    predMom <- pred.moments(x=x, mean=mean)
    result <- rep(NA_real_, length(theta))
    logwts <- log(support$weight)
    for (i in 1:length(theta))
      result[i] <- sum(exp(logwts + pnorm(q=theta[i], mean=predMom[,"mean"], sd=predMom[,"sd"], log.p=TRUE)))
    return(result)
  }

  qpredict <- function(p, x, mean=TRUE)
  {
    stopifnot(!missing(p), is.vector(p), is.numeric(p),
              all(is.finite(p)), all(p >= 0), all(p <= 1), 
              !missing(x), is.numeric(x), all(is.finite(x)),
              ((is.vector(x) && (length(x)==d))
               | (is.matrix(x) && ((ncol(x)==d) & ((nrow(x)==1) | (length(p)==1))))))
    quant <- function(pp, xx)
    {
      if (!is.finite(pp) || ((pp<0) | (pp>1))) {
        qu <- NA_real_
      } else if (pp==0) {
        qu <- -Inf
      } else if (pp==1) {
        qu <- Inf
      } else {
        mini <- -2
        while (ppredict(mini, xx, mean=mean) > pp) mini <- 2*mini
        maxi <- 2
        while (ppredict(maxi, xx, mean=mean) < pp) maxi <- 2*maxi
        qu <- uniroot(function(b){return(ppredict(b, xx, mean=mean) - pp)},
                      lower=mini, upper=maxi, tol=tol.uniroot)$root
      }
      return(qu)
    }
    if (is.matrix(x) && (nrow(x)>1)) {
      result <- apply(x, 1, function(xxx){quant(p,xxx)})
    } else {
      result <- apply(matrix(p,ncol=1), 1,
                      function(ppp){quant(ppp,as.vector(x))})
    }
    return(result)
  }

  rpredict <- function(n, x, mean=TRUE)
  {
    stopifnot(!missing(n), is.vector(n), is.numeric(n), length(n)==1,
              is.finite(n), n>=1, n==round(n), 
              !missing(x), is.vector(x), is.numeric(x),
              length(x) == d, all(is.finite(x)))
    predMom <- pred.moments(x=x, mean=mean)
    result <- NULL
    # determine numbers of samples from each mixture component:
    freq <- as.vector(stats::rmultinom(1, size=n, prob=support$weight))
    # draw actual samples component-wise:
    for (i in 1:length(support$weight)) {
      if (freq[i] > 0) {
        result <- c(result, rnorm(freq[i], mean=predMom[i,"mean"], sd=predMom[i,"sd"]))
      }
    }
    # re-shuffle rows:
    result <- result[sample(x=n, size=n, replace=FALSE)]
    return(result)
  }

  pred.interval <- function(level=0.95, x, mean=TRUE,
                            method=interval.type)
  {
    stopifnot(is.vector(level), is.numeric(level), length(level)==1,
              all(is.finite(level)), all(level > 0), all(level < 1), 
              !missing(x), is.numeric(x), all(is.finite(x)),
              ((is.vector(x) && (length(x)==d))
               | (is.matrix(x) && (ncol(x)==d))))
    if (level < 0.5) warning("\"level\" below 50% ... intended?")
      
    intfun <- function(xx)
    {
      if (method=="central") { # "central", equal-tailed interval:
        inter <- qpredict(p = c((1-level)/2, 1-(1-level)/2), xx, mean=mean)
      } else {                 # "shortest" interval:
        intwidth <- function(left)
        {
          pleft <- ppredict(left, xx, mean=mean)
          right <- qpredict(level + pleft, xx, mean=mean)
          return(right-left)
        }
        opti <- optimize(intwidth, lower=qpredict(p=(1-level)/50, xx, mean=mean), 
                         upper=qpredict(p=1-level, xx, mean=mean))$minimum
        inter <- c(opti, qpredict(p=level+ppredict(opti, xx, mean=mean), xx, mean=mean))
      }
      return(inter)
    }
    result <- t(apply(matrix(x, ncol=d), 1, intfun))
    if (nrow(result)==1) {
      result <- as.vector(result)
      names(result) <- c("lower","upper")
    } else {
      colnames(result) <- c("lower","upper")
    }
    attr(result, "interval.type") <- method
    return(result)
  }

  dshrink <- function(theta, which, log=FALSE)
  {
    stopifnot(!missing(theta), is.vector(theta), is.numeric(theta), all(is.finite(theta)),
              !missing(which), is.vector(which), length(which) == 1,
              (is.numeric(which) | is.character(which)))
    if (is.numeric(which)) {
      stopifnot(is.element(which, 1:k))
      idx <- which
    } else if (is.character(which)) {
      stopifnot(is.element(which, labels))
      idx <- which(which == labels)
    }
    shrinkMom <- shrink.moments(which=which)
    result <- rep(NA_real_, length(theta))
    logwts <- log(support$weight)
    for (i in 1:length(theta))
      result[i] <- sum(exp(logwts + dnorm(x=theta[i], mean=shrinkMom[,"mean"], sd=shrinkMom[,"sd"], log=TRUE)))
    if (log) result <- log(result)
    return(result)
  }

  pshrink <- function(theta, which)
  {
    stopifnot(!missing(theta), is.vector(theta), is.numeric(theta), all(is.finite(theta)),
              !missing(which), is.vector(which), length(which) == 1,
              (is.numeric(which) | is.character(which)))
    if (is.numeric(which)) {
      stopifnot(is.element(which, 1:k))
      idx <- which
    } else if (is.character(which)) {
      stopifnot(is.element(which, labels))
      idx <- which(which == labels)
    }
    shrinkMom <- shrink.moments(which=which)
    result <- rep(NA_real_, length(theta))
    logwts <- log(support$weight)
    for (i in 1:length(theta))
      result[i] <- sum(exp(logwts + pnorm(q=theta[i], mean=shrinkMom[,"mean"], sd=shrinkMom[,"sd"], log.p=TRUE)))
    return(result)
  }
  
  qshrink <- function(p, which)
  {
    stopifnot(!missing(p), is.vector(p), is.numeric(p),
              all(is.finite(p)), all(p >= 0), all(p <= 1), 
              !missing(which), is.vector(which), length(which) == 1,
              (is.numeric(which) | is.character(which)))
    if (is.numeric(which)) {
      stopifnot(is.element(which, 1:k))
    } else if (is.character(which)) {
      stopifnot(is.element(which, labels))
    }
    quant <- function(pp)
    {
      mini <- -2
      while (pshrink(mini, which) > pp) mini <- 2*mini
      maxi <- 2
      while (pshrink(maxi, which) < pp) maxi <- 2*maxi
      ur <- uniroot(function(s){return(pshrink(s, which) - pp)},
                    lower=mini, upper=maxi, tol=tol.uniroot)
      return(ur$root)      
    }
    result <- rep(NA,length(p))
    result[p==0] <- -Inf
    result[p==1] <- Inf
    proper <- ((p>0) & (p<1))
    if (any(proper)) result[proper] <- apply(matrix(p[proper],ncol=1), 1, quant)
    return(result)
  }

  rshrink <- function(n, which)
  {
    stopifnot(!missing(n), is.vector(n), is.numeric(n), length(n)==1,
              is.finite(n), n>=1, n==round(n), 
              !missing(which), is.vector(which), length(which) == 1,
              (is.numeric(which) | is.character(which)))
    if (is.numeric(which)) {
      stopifnot(is.element(which, 1:k))
      idx <- which
    } else if (is.character(which)) {
      stopifnot(is.element(which, labels))
      idx <- which(which == labels)
    }
    shrinkMom <- shrink.moments(which=which)
    result <- NULL
    # determine numbers of samples from each mixture component:
    freq <- as.vector(stats::rmultinom(1, size=n, prob=support$weight))
    # draw actual samples component-wise:
    for (i in 1:length(support$weight)) {
      if (freq[i] > 0) {
        result <- c(result, rnorm(freq[i], mean=shrinkMom[i,"mean"], sd=shrinkMom[i,"sd"]))
      }
    }
    # re-shuffle rows:
    result <- result[sample(x=n, size=n, replace=FALSE)]
    return(result)
  }
  
  shrink.interval <- function(level, which, method=interval.type)
  {
    stopifnot(!missing(level), is.vector(level), is.numeric(level), length(level)==1,
              all(is.finite(level)), all(level > 0), all(level < 1), 
              !missing(which), is.vector(which), length(which) == 1,
              (is.numeric(which) | is.character(which)))
    if (level < 0.5) warning("\"level\" below 50% ... intended?")
    if (method=="central") { # "central", equal-tailed interval:
      result <- qshrink(p = c((1-level)/2, 1-(1-level)/2), which)
    } else {                 # "shortest" interval:
      intwidth <- function(left)
      {
        pleft <- pshrink(left, which)
        right <- qshrink(level + pleft, which)
        return(right-left)
      }
      opti <- optimize(intwidth, lower=qshrink(p=(1-level)/50, which), 
                       upper=qshrink(p=1-level, which))$minimum
      result <- c(opti, qshrink(p=level+pshrink(opti, which), which))
    }
    attr(result, "interval.type") <- method
    return(result)
  }
  
  marglikeli <- function(tau, beta, which.beta)
  {
    if (!missing(tau)) {
      stopifnot(is.numeric(tau), length(tau)==1)
    } else {
      tau <- NA_real_
    }
    if (!missing(beta)) {
      betaMatrix <- beta.convert(beta, which.beta, d, betanames)
      stopifnot(nrow(betaMatrix)==1)
      betaVector <- as.vector(betaMatrix)
    } else {
      betaVector <- rep(NA_real_, d)
    }
    
    marglikfun <- function(tau, beta)
    # compute CONDITIONAL (marginal) likelihood (for given tau) 
    {
      if (tau<0) {
        margl <- 0.0
      } else {
        # (marginal) prior mean and cov. of beta:
        betaMean <- betaprior$mean
        betaCov  <- betaprior$covariance
        idx <- is.finite(beta) # (the conditioned-upon beta elements)
        if (any(idx)) { # determine conditional prior moments
          if (!all(idx)) {
            Sigma12 <- betaCov[!idx,idx]
            Sigma22 <- betaCov[idx,idx]
            Sigma22inv <- solve(betaCov[idx,idx])
            betaMean[!idx] <- betaMean[!idx] + Sigma12 %*% Sigma22inv  %*% (beta[idx]-betaMean[idx])
            betaCov[!idx,!idx] <- betaCov[!idx,!idx] - Sigma12 %*% Sigma22inv %*% t(Sigma12)
          }
          betaMean[idx] <- beta[idx]
          betaCov[idx,] <- 0.0
          betaCov[,idx] <- 0.0
        }
        if (all(beta.prior.proper | idx)) {
          ppMean <- X %*% betaMean
          ppCov  <- diag(sigma^2 + tau^2) + X %*% betaCov %*% t(X)
          margl <- mvtnorm::dmvnorm(y, mean=ppMean, sigma=ppCov)
        } else {
          margl <- 0.0
        }
      }
      return(margl)
    }

    if (is.finite(tau)) { # conditional on tau
      result <- marglikfun(tau, betaVector)
    } else {              # marginalized over tau
      if (tau.prior.proper) {
        integrand <- function(t)
        {
          return(apply(matrix(t,ncol=1), 1,
                       function(tt){marglikfun(tt,betaVector) * dprior(tau=tt)}))
        }
        int <- integrate(integrand, lower=0, upper=Inf,
                         rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate,
                         stop.on.error=FALSE)
        if ((int$message == "OK") && (int$value > 0))
          result <- int$value
        else
          result <- 0.0
      } else {
        result <- 0.0
      }
    }    
    return(result)
  }
  
  conditionalmoments <- function(tau)
  # posterior moments (mean and covariance) of regression parameters (beta)
  # conditional on a given heterogeneity (tau) value.
  {
    stopifnot(length(tau)==1, is.finite(tau), tau>=0)
    if (!is.null(Sigmap)) { # set up block diagonal matrix:
      sigmaTauInv <- cbind(rbind(diag(1 / (sigma^2 + tau^2), nrow=k, ncol=k),
                                 matrix(0.0, nrow=length(yp), ncol=k)),
                           rbind(matrix(0.0, nrow=k, ncol=length(yp)),
                                 solve(Sigmap)))
    } else {
      sigmaTauInv <- diag(1 / (sigma^2 + tau^2), nrow=k, ncol=k)
    }
    Vbeta       <- solve(t(rbind(X,Xp)) %*% sigmaTauInv %*% rbind(X,Xp))
    betaHat     <- as.vector(Vbeta %*% t(rbind(X,Xp)) %*% sigmaTauInv %*% c(y,yp))
    return(list("mean"=betaHat, "covariance"=Vbeta))
  }

  pred.moments <- function(tau, x, mean=TRUE)
  # posterior predictive moments (mean and std. dev.)
  # conditional on given heterogeneity ("tau") values
  # and covariable vector "x".
  {
    stopifnot(!missing(x), is.vector(x), is.numeric(x),
              length(x) == d, all(is.finite(x)))
    if (missing(tau)) { # return moments based on "support" object
      tau        <- support$tau
      meanMatrix <- support$mean
      covArray   <- support$covariance
    } else {            # return moments based on supplied "tau" argument
      stopifnot(is.vector(tau), all(is.finite(tau)), all(tau>=0))
      meanMatrix <- matrix(NA_real_, nrow=length(tau), ncol=d,
                           dimnames=list(NULL, betanames))
      covArray   <- array(NA_real_, dim=c(length(tau), d, d),
                          dimnames=list(NULL, betanames, betanames))
      for (i in 1:length(tau)) {
        cm <- conditionalmoments(tau=tau[i])
        meanMatrix[i,] <- cm$mean
        covArray[i,,]  <- cm$covariance
      }
    }
    # derive predictive moments:
    meanvec <- as.vector(meanMatrix %*% x)
    sdvec   <- rep(NA_real_, nrow(meanMatrix))
    for (i in 1:nrow(meanMatrix))
      sdvec[i] <- as.vector(t(x) %*% covArray[i,,] %*% x)
    if (!mean) sdvec <- sdvec + tau^2
    sdvec <- sqrt(sdvec)
    return(cbind("mean"=meanvec, "sd"=sdvec))
  }

  shrink.moments <- function(tau, which)
  # shrinkage moments (mean and std. dev.)
  # conditional on given heterogeneity ("tau") values.
  {
    stopifnot(!missing(which), is.vector(which), length(which) == 1,
              (is.numeric(which) | is.character(which)))
    # derive (numerical) study identifier "idx":
    if (is.numeric(which)) {
      stopifnot(is.element(which, 1:k))
      idx <- which
    } else if (is.character(which)) {
      stopifnot(is.element(which, labels))
      idx <- which(which == labels)
    }
    if (missing(tau)) { # return moments based on "support" object
      tauVec  <- support$tau
      predMom <- pred.moments(x=X[idx,])
    } else {            # return moments based on supplied "tau" argument
      stopifnot(is.vector(tau), all(is.finite(tau)), all(tau>=0))
      tauVec  <- tau
      predMom <- pred.moments(tau=tau, x=X[idx,])
    }
    # derive shrinkage moments;
    # "inverse variance" weights:
    ivw1 <- sigma[idx]^-2 / (sigma[idx]^-2 + tauVec^-2)
    ivw2 <- tauVec^-2 / (sigma[idx]^-2 + tauVec^-2)
    ivw2[tauVec==0] <- 1.0
    # shrinkage means:
    meanvec <- ivw1*y[idx] + ivw2*predMom[,"mean"]
    # shrinkage stdevs:
    sdvec   <- sqrt((1 / (sigma[idx]^-2 + tauVec^-2)) + (ivw2*predMom[,"sd"])^2)
    return(cbind("mean"=meanvec, "sd"=sdvec))
  }

  discretize <- function(delta=0.01, epsilon=0.0001)
  {
    divergence <- function(tau1, tau2)
    {
      cm1 <- conditionalmoments(tau=tau1)
      cm2 <- conditionalmoments(tau=tau2)
      return(kldiv(mu1=cm1$mean, mu2=cm2$mean,
                   sigma1=cm1$covariance, sigma2=cm2$covariance,
                   symmetrized=TRUE))
    }
    tau <- 0.0  # (first bin's reference point)
    maxtau <- qposterior(tau.p = 1 - min(c(epsilon, 1-epsilon))) # (end of relevant range)
    # search for upper bin margin:
    upper  <- maxtau
    diverg <- divergence(tau, upper)
    if (diverg > delta) {
      ur <- uniroot(function(t){divergence(tau,t)-delta},
                    lower=tau, upper=upper,
                    f.lower=-delta, f.upper=diverg-delta,
                    tol=tol.uniroot)
      tau <- ur$root
    } else {
      tau <- upper
    }
    prob1 <- 0.0
    prob2 <- pposterior(tau=tau)
    # store result for 1st bin:
    result <- matrix(c(0.0, prob2-prob1),
                     nrow=1, ncol=2,
                     dimnames=list(NULL,c("tau","weight")))
    # determine following bins (2,...):
    bin <- 2
    while ((tau < maxtau) | (bin <= 2)) {  # (enforce at least 2 support points)
      laststep <- tau - result[nrow(result),"tau"]
      result <- rbind(result, rep(0,2))
      # determine new bin's reference point:
      upper <- tau + 2*laststep
      diverg <- divergence(tau, upper)
      while ((diverg <= delta) & (upper<=maxtau)) {
        upper  <- upper + laststep
        diverg <- divergence(tau, upper)
      }
      if (diverg <= delta) { # early stop; set last reference & margin "manually":
        tau <- tau + max(c(laststep, (maxtau-tau) / 2))
        result[bin,"tau"] <- tau
        tau <- max(c(tau + laststep), maxtau)
        # determine bin's weight:
        prob1 <- prob2
        prob2 <- pposterior(tau=tau)
        result[bin,"weight"] <- max(c(0.0, prob2-prob1))
      } else {               # search for reference point & margin:
        ur <- uniroot(function(t){return(divergence(tau, t)-delta)},
                      lower=tau, upper=upper,
                      f.lower=-delta, f.upper=diverg-delta,
                      tol=tol.uniroot)
        laststep <- ur$root - tau
        tau <- ur$root
        result[bin,"tau"] <- tau
        # determine bin's upper bound:
        upper <- tau + 2*laststep
        diverg <- divergence(tau, upper)
        while ((diverg <= delta) & (upper <= maxtau)) {
          upper <- upper + laststep
          diverg <- divergence(tau, upper)
        }
        if (diverg <= delta) { # early stop; set last margin "manually":
          tau <- max(c(tau + laststep), maxtau)
          # determine bin's weight:
          prob1 <- prob2
          prob2 <- pposterior(tau=tau)
          result[bin,"weight"] <- max(c(0.0, prob2-prob1))
        } else {               # search for reference point & margin:
          ur <- uniroot(function(t){return(divergence(tau, t)-delta)},
                        lower=tau, upper=upper,
                        f.lower=-delta, f.upper=diverg-delta,
                        tol=tol.uniroot)
          tau <- ur$root
          # determine bin's weight:
          prob1 <- prob2
          prob2 <- pposterior(tau=tau)
          result[bin,"weight"] <- max(c(0.0, prob2-prob1))
        }
      }
      # sanity check (to catch possible uniroot() failure):
      if (result[bin,"tau"] <= result[bin-1,"tau"]) {
        warning("DIRECT grid setup appears to have failed.")
        tau <- maxtau + 1.0
      }
      bin <- bin+1
    }
    # re-normalize weights (if necessary):
    sw <- sum(result[,"weight"])
    if (sw <= 0) {
      warning("DIRECT grid setup returned strange weights.")
      result[,"weight"] <- rep(1/nrow(result), nrow(result))
      sw <- 1.0
    }
    if (sw != 1.0)
      result[,"weight"] <- result[,"weight"] / sw
    return(result)
  }

  ############################################################
  # derive set of support points via DIRECT method:
  direct <- discretize(delta=delta, epsilon=epsilon)
  support <- list("weight" = direct[,"weight"],
                  "tau"    = direct[,"tau"],
                  "mean"   = NULL, "covariance"=NULL)
  support$mean       <- matrix(NA_real_, nrow=length(support$weight), ncol=d,
                               dimnames=list(NULL, betanames))
  support$covariance <- array(NA_real_, dim=c(length(support$weight), d, d),
                              dimnames=list(NULL, betanames, betanames))
  for (i in 1:length(support$weight)) {
    cm <- conditionalmoments(tau=support$tau[i])
    support$mean[i,]        <- cm$mean
    support$covariance[i,,] <- cm$covariance
  }
  rm(list=c("direct", "cm"))

  # derive posterior summary stats:
  sumstats <- matrix(NA_real_, nrow=6, ncol=1+d,
                     dimnames=list(c("mode", "median", "mean", "sd",
                                     "95% lower", "95% upper"),
                                   c("tau", betanames)))
  # heterogeneity (tau) moments:
  expectation <- try(integrate(function(x){return(dposterior(tau=x)*x)},
                               lower=0, upper=Inf,
                               rel.tol=rel.tol.integrate,
                               abs.tol=abs.tol.integrate)$value, silent=TRUE)
  if (inherits(expectation,"try-error")) {
    expectation <- NA
    variance <- NA
  } else {
    variance <- try(integrate(function(x){return(dposterior(tau=x)*(x-expectation)^2)},
                              lower=0, upper=Inf,
                              rel.tol=rel.tol.integrate,
                              abs.tol=abs.tol.integrate)$value, silent=TRUE)
    if (inherits(variance, "try-error"))
      variance <- NA
  }
  sumstats[c("mean","sd"), "tau"] <- c(expectation, sqrt(variance))
  rm(list=c("expectation", "variance"))
  # heterogeneity (tau) median:
  sumstats["median","tau"] <- qposterior(tau.p=0.5)
  # regression parameters' (beta) (marginal) medians:
  for (i in 1:d)
    sumstats["median",betanames[i]] <- qposterior(beta.p=0.5, which=i)
  # shortest credible intervals:
  sumstats[c("95% lower","95% upper"),"tau"] <- post.interval(tau.level=0.95,
                                                              method=interval.type)
  for (i in 1:d)
    sumstats[c("95% lower","95% upper"),betanames[i]] <- post.interval(beta.level=0.95,
                                                                       which.beta=i,
                                                                       method=interval.type)
  # heterogeneity (tau) mode:
  if (! is.finite(dposterior(tau=0))) {
    sumstats["mode","tau"] <- 0
  } else {
    maxi <- optimize(function(x){return(dposterior(tau=x))},
                     lower=0, upper=qposterior(tau.p=0.9), maximum=TRUE)
    if ((maxi$maximum <= 2 * .Machine$double.eps^0.25)
        && (dposterior(tau=0) > maxi$objective)) maxi <- list("maximum"=0)
    sumstats["mode","tau"] <- maxi$maximum
    rm(list="maxi")
  }

  # regression parameters' (beta) moments
  betamean <- as.vector(t(support$mean) %*% support$weight)
  # (co-)variances:
  betacov <- matrix(NA_real_, nrow=d, ncol=d)
  for (i in 1:d) {
    for (j in i:d) {
      # compute weighted averages:
      betacov[i,j] <- t(support$covariance[,i,j]) %*% support$weight
      if (j==i) { # main diagonal (variances):
        betacov[i,i] <- betacov[i,i] + t((support$mean[,i]-betamean[i])^2) %*% support$weight
      } else {    # off-diagonal (covariances):
        betacov[i,j] <- betacov[i,j] + t((support$mean[,i]-betamean[i])*(support$mean[,j]-betamean[j])) %*% support$weight
        betacov[j,i] <- betacov[i,j]
      }
    }
  }
  sumstats["mean",betanames] <- betamean
  sumstats["sd",  betanames] <- sqrt(diag(betacov))
  beta.moments <- list("mean"=betamean, "covariance"=betacov)
  names(beta.moments$mean) <- rownames(beta.moments$covariance) <- colnames(beta.moments$covariance) <- betanames
  rm(list=c("betamean", "betacov"))
  
  # regression parameters' (beta) (marginal) modes:
  for (i in 1:d) {
    maxi <- optimize(function(b){dposterior(beta=b, which.beta=i)},
                     lower=sumstats["mean",betanames[i]] - 3*sumstats["sd",betanames[i]],
                     upper=sumstats["mean",betanames[i]] + 3*sumstats["sd",betanames[i]],
                     maximum=TRUE)
    sumstats["mode",betanames[i]] <- maxi$maximum
  }
  rm(list="maxi")

  # MAP estimates:
  map.estimate <- matrix(NA_real_, ncol=d+1, nrow=2,
                         dimnames=list(c("joint","marginal"),
                                       c("tau", betanames)))
  map.estimate["marginal",] <- sumstats["mode",]
  # global optimization:
  opti <- optim(sumstats["median",],
                function(x){return(-dposterior(tau=x[1], beta=x[2:(d+1)], which.beta=1:d, log=TRUE))})
  # optimization with tau=0 (to potentially find mode at margin):
  if (is.finite(dposterior(tau=0.0, beta=sumstats["median",-1], which.beta=1:d, log=TRUE))){
    if (d>1)
      opti0 <- optim(sumstats["median",-1],
                     function(x){return(-dposterior(tau=0.0, beta=x[1:d], which.beta=1:d, log=TRUE))})
    else
      opti0 <- optim(sumstats["median",-1],
                     function(x){return(-dposterior(tau=0.0, beta=x[1:d], which.beta=1:d, log=TRUE))},
                     method="Brent", lower=sumstats["95% lower",2], upper=sumstats["95% upper",2])
  } else {
    opti0 <- list("value" = Inf)
  }
  # return the better of the two optima:
  if (opti0$value < opti$value)
    map.estimate["joint",] <- c(0.0, opti0$par)
  else 
    map.estimate["joint",] <- opti$par
  
  # compute "shrinkage" estimates of theta[i]:
  shrink <- matrix(NA, nrow=8, ncol=k,
                   dimnames=list(c("y","sigma","mode", "median", "mean","sd", "95% lower", "95% upper"),
                                 labels))
  shrink["y",] <- y
  shrink["sigma",] <- sigma
  for (i in 1:k) {
    shrink["median",i] <- qshrink(p=0.5, which=i)
    shrink[c("95% lower", "95% upper"),i] <- shrink.interval(level=0.95, which=i)
    maxi <- optimize(function(x){return(dshrink(theta=x, which=i))},
                     lower=shrink["95% lower",i],
                     upper=shrink["95% upper",i], maximum=TRUE)
    shrink["mode",i] <- maxi$maximum
    # compute moments; predictive moments:
    meanvec <- as.vector(support$mean %*% X[i,])
    sdvec   <- rep(NA_real_, length(support$weight))
    for (j in 1:length(support$weight))
      sdvec[j] <- as.vector(t(X[i,]) %*% support$covariance[j,,] %*% X[i,])
    sdvec <- sqrt(sdvec)
    # derive shrinkage moments:
    # "inverse variance" weights:
    ivw1 <- sigma[i]^-2 / (sigma[i]^-2 + support$tau^-2)
    ivw2 <- support$tau^-2 / (sigma[i]^-2 + support$tau^-2)
    ivw2[support$tau==0] <- 1.0
    # (conditional) shrinkage means:
    meanvec <- ivw1*y[i] + ivw2*meanvec
    # (conditional) shrinkage stdevs.:
    sdvec   <- sqrt((1 / (sigma[i]^-2 + support$tau^-2)) + (ivw2*sdvec)^2)
    #print(cbind("mean"=meanvec, "sd"=sdvec))
    shrink["mean",i] <- sum(meanvec * support$weight)
    shrink["sd",i]   <- sqrt(sum(((meanvec-shrink["mean",i])^2 + sdvec^2) * support$weight))
  }

  # compute ML estimate(s); joint:
  ml.estimate <- matrix(NA_real_, nrow=2, ncol=d+1,
                        dimnames=list(c("joint","marginal"), c("tau",betanames)))
  opti <- optim(sumstats["median",],
                function(x){return(-marglikeli(tau=x[1], beta=x[2:(d+1)], which.beta=1:d))})
  if (d>1)
    opti0 <- optim(sumstats["median",-1],
                   function(x){return(-marglikeli(tau=0.0, beta=x[1:d], which.beta=1:d))})
  else
    opti0 <- optim(sumstats["median",-1],
                   function(x){return(-marglikeli(tau=0.0, beta=x[1:d], which.beta=1:d))},
                   method="Brent", lower=sumstats["95% lower",2], upper=sumstats["95% upper",2])
  if (opti0$value < opti$value) {
    ml.estimate["joint",] <- c(0.0, opti0$par)
    ml.value              <- -opti0$value
  } else {
    ml.estimate["joint",] <- opti$par
    ml.value              <- -opti$value
  }
  # compute marginal ML estimates:
  if (all(beta.prior.proper)) {
    upper <- 2 * sumstats["median","tau"]
    while (marglikeli(tau=2*upper) > marglikeli(tau=upper)) upper <- upper * 2
    maxi <- optimize(function(x){return(marglikeli(tau=x))},
                     lower=0, upper=upper*2, maximum=TRUE)
    ml.estimate["marginal","tau"] <- maxi$maximum
  }
  if (tau.prior.proper) {
    for (i in 1:d) {
      if ((d==1) || all(beta.prior.proper[-i])) {
        upper <- ml.estimate["joint",i+1] + 1
        while (marglikeli(beta=upper+1, which.beta=i) > marglikeli(beta=upper, which.beta=i))
          upper <- upper + 1
        lower <- ml.estimate["joint",i+1] - 1
        while (marglikeli(beta=lower-1, which.beta=i) > marglikeli(beta=lower, which.beta=i))
          lower <- lower - 1
        maxi <- optimize(function(x){return(marglikeli(beta=x, which.beta=i))},
                         lower=lower-1, upper=upper+1, maximum=TRUE)
        ml.estimate["marginal",i+1] <- maxi$maximum
      }
    }
  }
  
  # compute marginal likelihood:
  if (tau.prior.proper & all(beta.prior.proper))
    marglik <- marglikeli()
  else {
    marglik <- NA_real_
    if (!tau.prior.proper & !all(beta.prior.proper))
      attr(marglik, "NA.reason") <- "improper heterogeneity and effect priors"
    else if (!tau.prior.proper)
      attr(marglik, "NA.reason") <- "improper heterogeneity prior"
    else
      attr(marglik, "NA.reason") <- "improper effect prior"
  }
  
  # compute Bayes factors:
  bayesfactor <- matrix(NA_real_, nrow=2, ncol=d+1,
                        dimnames=list(c("actual", "minimum"), c("tau=0", paste0(betanames,"=0"))))
  # "actual" Bayes factors:
  if (is.finite(marglik)) {
    bayesfactor["actual","tau=0"]  <- marglikeli(tau=0) / marglik
    for (i in 1:d) {
      bayesfactor["actual",i+1]  <- marglikeli(beta=0, which.beta=i) / marglik
    }
  }
  # "minimum" Bayes factors:
  if (all(beta.prior.proper))
    bayesfactor["minimum","tau=0"]  <- marglikeli(tau=0) / marglikeli(tau=ml.estimate["marginal","tau"])
  if (tau.prior.proper) {
    for (i in 1:d) {
      if ((d==1) || all(beta.prior.proper[-i])) {
        bayesfactor["minimum",i+1]  <- marglikeli(beta=0, which.beta=i) / marglikeli(beta=ml.estimate["marginal",i+1], which.beta=i)
      }
    }
  }
  
  # wrap up:
  ptm <- proc.time()[1] - ptm[1]  
  result <- list("y"=y, "sigma"=sigma, "X"=X,
                 "k"=k, "d"=d,
                 "labels"=labels, "variables"=betanames,
                 "tau.prior"           = tau.prior,
                 "tau.prior.proper"    = tau.prior.proper,
                 "beta.prior"          = betaprior,
                 "beta.prior.proper"   = beta.prior.proper,
                 "dprior"              = dprior,
                 "likelihood"          = marglikeli,
                 "dposterior"          = dposterior,
                 "pposterior"          = pposterior,
                 "qposterior"          = qposterior,
                 "rposterior"          = rposterior,
                 "post.interval"       = post.interval,
                 "dpredict"            = dpredict,
                 "ppredict"            = ppredict,
                 "qpredict"            = qpredict,
                 "rpredict"            = rpredict,
                 "pred.interval"       = pred.interval,
                 "dshrink"             = dshrink,
                 "pshrink"             = pshrink,
                 "qshrink"             = qshrink,
                 "rshrink"             = rshrink,
                 "shrink.interval"     = shrink.interval,
                 "post.moments"        = conditionalmoments,
                 "pred.moments"        = pred.moments,
                 "shrink.moments"      = shrink.moments,
                 "summary"             = sumstats,
                 "interval.type"       = interval.type,
                 "ML"                  = ml.estimate,
                 "MAP"                 = map.estimate,
                 "theta"               = shrink,
                 "marginal.likelihood" = marglik,
                 "bayesfactor"         = bayesfactor,
                 "beta.moments"        = beta.moments,
                 "support"             = support,
                 "delta"               = delta,
                 "epsilon"             = epsilon,
                 "rel.tol.integrate"   = rel.tol.integrate,
                 "abs.tol.integrate"   = abs.tol.integrate,
                 "tol.uniroot"         = tol.uniroot,
                 "call"                = match.call(expand.dots=FALSE),
                 "init.time"           = c("seconds"=unname(ptm)))
  class(result) <- "bmr"
  return(result)  
}


########################################


bmr.escalc <- function(y, labels=NULL, ...)
# apply "bmr()" to an "escalc" object
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
  result <- bmr.default(y=as.vector(y[,var.names[1]]),
                        sigma=sqrt(as.vector(y[,var.names[2]])),
                        labels=labels, ...)
  result$call <- match.call(expand.dots=FALSE)  
  return(result)
}


########################################


print.bmr <- function(x,...)
# print a short summary
{
  cat(" 'bmr' object.\n")
  cat(paste("\n",x$k," estimates:\n", sep=""))
  if (length(x$labels)>10)
    cat(paste(c(x$labels[1:10], "..."), collapse=", "))
  else
    cat(paste(x$labels[1:length(x$labels)], collapse=", "))
  cat("\n")
  if (x$d == 1) cat(paste0("\n",x$d," regression parameter:\n"))
  else cat(paste0("\n",x$d," regression parameters:\n"))
  cat(paste(x$variables, collapse=", "))
  cat("\n")
  cat(paste0("\ntau prior (", ifelse(x$tau.prior.proper, "proper", "improper"), "):\n"))
  if (is.element("bayesmeta.label", names(attributes(x$tau.prior))))
    cat(paste(attributes(x$tau.prior)[["bayesmeta.label"]], "\n"))
  else
    print(x$tau.prior)
  if (any(x$beta.prior.proper)) {
    cat(paste0("\nbeta prior (",
               ifelse(all(x$beta.prior.proper), "proper", "partly proper"), "):\n"))
    print(x$beta.prior)
  } else {
    cat("\nbeta prior: (improper) uniform\n")
  }
  cat("\nMAP estimates:\n")
  print(x$MAP)
  cat("\nmarginal posterior summary:\n")
  print(x$summary)
  if (x$interval.type=="shortest")
    cat("\n(quoted intervals are shortest credible intervals.)\n")
  else
    cat("\n(quoted intervals are central, equal-tailed credible intervals.)\n")      
  invisible(x)
}


summary.bmr <- function(object, X.mean, X.prediction, ...)
# compute "summary.bmr" object
{
  ############################################################
  # default treatment of "X.mean" argument:
  if (missing(X.mean) || (is.null(X.mean) || all(is.na(X.mean)))) {
      X.mean <- NULL
      meanNumber <- 0
  } else {
    stopifnot(is.numeric(X.mean),
              (is.matrix(X.mean) && (ncol(X.mean) == object$d))
              || (is.vector(X.mean) && (length(X.mean) == object$d)))
    if (is.vector(X.mean)) {
      X.mean <- matrix(X.mean, nrow=1, ncol=object$d)
    }
    if (is.null(colnames(X.mean))) colnames(X.mean) <- object$variables      
    meanNumber <- nrow(X.mean)
  }
  # default treatment of "X.prediction" argument:
  if (missing(X.prediction) || (is.null(X.prediction) || all(is.na(X.prediction)))) {
      X.prediction <- NULL
      predNumber <- 0
  } else {
    stopifnot(is.numeric(X.prediction),
              (is.matrix(X.prediction) && (ncol(X.prediction) == object$d))
              | (is.vector(X.prediction) && (length(X.prediction) == object$d)))
    if (is.vector(X.prediction)) {
      X.mean <- matrix(X.prediction, nrow=1, ncol=object$d)
    }
    if (is.null(colnames(X.prediction))) colnames(X.prediction) <- object$variables      
    predNumber <- nrow(X.prediction)
  }
  ############################################################
  # compute means and/or predictions:
  summary.mean <- summary.pred <- NULL
  if (meanNumber > 0) { # compute mean estimates & CIs:
    summary.mean <- matrix(NA_real_, nrow=meanNumber, ncol=6,
                           dimnames=list(rownames(X.mean),
                                         c("mode", "median", "mean", "sd",
                                           "95% lower", "95% upper")))
    summary.mean[,"median"] <- object$qpredict(p=0.5, x=X.mean, mean=TRUE)
    summary.mean[,c("95% lower", "95% upper")] <- object$pred.interval(x=X.mean, mean=TRUE,
                                                                       method=object$interval.type)
    # compute marginal moments & modes:
    for (i in 1:meanNumber) {
      # first compute expectation:
      expect <- integrate(function(b){return(object$dpredict(theta=b, x=X.mean[i,], mean=TRUE) * b)},
                          lower=-Inf, upper=Inf,
                          rel.tol=object$rel.tol.integrate, abs.tol=object$abs.tol.integrate,
                          stop.on.error=FALSE)
      if (expect$message == "OK") {
        # if expectation is finite, store and proceed with variance:
        summary.mean[i,"mean"] <- expect$value
        vari <- integrate(function(b){return(object$dpredict(theta=b, x=X.mean[i,], mean=TRUE) * (b-expect$value)^2)},
                          lower=-Inf, upper=Inf,
                          rel.tol=object$rel.tol.integrate, abs.tol=object$abs.tol.integrate,
                          stop.on.error=FALSE)
        if (vari$message == "OK") {
          summary.mean[i,"sd"] <- sqrt(vari$value)
        }
      }
      maxi <- optimize(function(b){return(object$dpredict(theta=b, x=X.mean[i,], mean=TRUE))},
                       lower=summary.mean[,"95% lower"], upper=summary.mean[,"95% upper"],
                       maximum=TRUE)
      summary.mean[i,"mode"] <- maxi$maximum
    }
  }
  if (predNumber > 0) { # compute predictions & CIs:
    summary.pred <- matrix(NA_real_, nrow=predNumber, ncol=6,
                           dimnames=list(rownames(X.prediction),
                                         c("mode", "median", "mean", "sd",
                                           "95% lower", "95% upper")))
    summary.pred[,"median"] <- object$qpredict(p=0.5, x=X.prediction, mean=FALSE)
    summary.pred[,c("95% lower", "95% upper")] <- object$pred.interval(x=X.prediction, mean=FALSE,
                                                                       method=object$interval.type)
    # compute marginal moments & modes:
    for (i in 1:predNumber) {
      # first compute expectation:
      expect <- integrate(function(b){return(object$dpredict(theta=b, x=X.prediction[i,], mean=FALSE) * b)},
                          lower=-Inf, upper=Inf,
                          rel.tol=object$rel.tol.integrate, abs.tol=object$abs.tol.integrate,
                          stop.on.error=FALSE)
      if (expect$message == "OK") {
        # if expectation is finite, store and proceed with variance:
        summary.pred[i,"mean"] <- expect$value
        vari <- integrate(function(b){return(object$dpredict(theta=b, x=X.prediction[i,], mean=FALSE) * (b-expect$value)^2)},
                          lower=-Inf, upper=Inf,
                          rel.tol=object$rel.tol.integrate, abs.tol=object$abs.tol.integrate,
                          stop.on.error=FALSE)
        if (vari$message == "OK") {
          summary.pred[i,"sd"] <- sqrt(vari$value)
        }
      }
      maxi <- optimize(function(b){return(object$dpredict(theta=b, x=X.prediction[i,], mean=FALSE))},
                       lower=summary.pred[,"95% lower"], upper=summary.pred[,"95% upper"],
                       maximum=TRUE)
      summary.pred[i,"mode"] <- maxi$maximum
    }
  }
  result <- list("bmr"          = object,
                 "call"         = match.call(expand.dots=FALSE),
                 "X.mean"       = X.mean,
                 "X.prediction" = X.prediction,
                 "mean"         = summary.mean,
                 "prediction"   = summary.pred)
  class(result) <- "summary.bmr"
  return(result)
}


print.summary.bmr <- function(x, ...)
# print a longer summary
{
  cat(" 'bmr' object.\n")
  cat(paste("\n",x$bmr$k," estimates:\n", sep=""))
  if (length(x$bmr$labels)>10)
    cat(paste(c(x$bmr$labels[1:10], "..."), collapse=", "))
  else
    cat(paste(x$bmr$labels[1:length(x$bmr$labels)], collapse=", "))
  cat("\n")
  if (x$bmr$d == 1) cat(paste0("\n",x$bmr$d," regression parameter:\n"))
  else cat(paste0("\n",x$bmr$d," regression parameters:\n"))
  cat(paste(x$bmr$variables, collapse=", "))
  cat("\n")
  cat(paste0("\ntau prior (", ifelse(x$bmr$tau.prior.proper, "proper", "improper"), "):\n"))
  if (is.element("bayesmeta.label", names(attributes(x$bmr$tau.prior))))
    cat(paste(attributes(x$bmr$tau.prior)[["bayesmeta.label"]], "\n"))
  else
    print(x$bmr$tau.prior)
  if (any(x$bmr$beta.prior.proper)) {
    cat(paste0("\nbeta prior (",
               ifelse(all(x$bmr$beta.prior.proper), "proper", "partly proper"), "):\n"))
    print(x$bmr$beta.prior)
  } else {
    cat("\nbeta prior: (improper) uniform\n")
  }
  cat("\nMAP estimates:\n")
  print(x$bmr$MAP)
  cat("\nmarginal posterior summary:\n")
  print(x$bmr$summary)
  ############################################################
  # print means & predictions:
  if (!is.null(x$mean)) {
    cat("\nmean estimates:\n")
    print(cbind(x$X.mean,
                x$mean[,c("median", "95% lower", "95% upper"), drop=FALSE]))
  }
  if (!is.null(x$prediction)) {
    cat("\npredictions:\n")
    print(cbind(x$X.prediction,
                x$prediction[,c("median", "95% lower", "95% upper"), drop=FALSE]))
  }
  ############################################################
  if (x$bmr$interval.type=="shortest")
    cat("\n(quoted intervals are shortest credible intervals.)\n")
  else
    cat("\n(quoted intervals are central, equal-tailed credible intervals.)\n")      
  if (any(is.finite(x$bmr$bayesfactor))) {
    cat("\nBayes factors:\n")
    print(x$bmr$bayesfactor)
  }
  invisible(x)
}


########################################


plot.bmr <- function(x, nrow, ncol, prior=FALSE, ...)
# draw posterior densities for each parameter
{
  par.mfrow <- par("mfrow")
  on.exit(par(mfrow=par.mfrow))
  # determine plot layout:
  if (missing(nrow) & !missing(ncol)) {
    stopifnot(length(ncol)==1, ncol>=1, ncol==round(ncol))
    nrow <- ceiling((x$d+1) / ncol)
  } else if (!missing(nrow) & missing(ncol)) {
    stopifnot(length(nrow)==1, nrow>=1, nrow==round(nrow))
    ncol <- ceiling((x$d+1) / nrow)      
  } else if (missing(nrow) & missing(ncol)) {
    nrow <- floor(sqrt(x$d+1))
    ncol <- ceiling((x$d+1) / nrow)      
  } else {
    stopifnot(length(nrow)==1, nrow>=1, nrow==round(nrow),
              length(ncol)==1, ncol>=1, ncol==round(ncol),
              nrow * ncol >= (x$d+1))
  }
  par(mfrow=c(nrow,ncol))
  xrange <- c(0, x$qposterior(tau=0.9990))
  tau <- seq(xrange[1], xrange[2], le=400)
  dens <- x$dposterior(tau=tau)
  maxdens <- max(dens[is.finite(dens)])
  dens[!is.finite(dens)] <- 2 * maxdens
  plot(xrange + c(0,-1)*diff(xrange)/20, c(0, maxdens), type="n",
       main="heterogeneity", xlab=expression(tau), ylab="density")
  # draw left tail:
  if (x$summary["95% lower","tau"] > 0) {
    idx <- which(tau < x$summary["95% lower","tau"])
    polygon(c(0, tau[idx], rep(x$summary["95% lower","tau"],2)),
            c(0, dens[idx], x$dposterior(tau=x$summary["95% lower","tau"]), 0),
            col="grey90", border=NA)
  }
  # draw central bulk:
  idx <- which((tau > x$summary["95% lower","tau"])
               & (tau < x$summary["95% upper","tau"]))
  poly.tau  <- c(rep(x$summary["95% lower","tau"],2), tau[idx], rep(x$summary["95% upper","tau"],2))
  poly.dens <- c(0, x$dposterior(tau=x$summary["95% lower","tau"]),
                 dens[idx], x$dposterior(tau=x$summary["95% upper","tau"]), 0)
  poly.dens[!is.finite(poly.dens)] <- maxdens
  polygon(poly.tau, poly.dens, col="grey80", border=NA)
  # draw right tail:
  idx <- which(tau > x$summary["95% upper","tau"])
  polygon(c(rep(x$summary["95% upper","tau"],2), tau[idx], rep(max(tau),2)),
          c(0, x$dposterior(tau=x$summary["95% upper","tau"]),
            dens[idx], dens[length(tau)], 0),
          col="grey90", border=NA)
  # draw axes:
  abline(h=0, v=0, col="grey40")
  # draw median:
  lines(rep(x$summary["median","tau"],2),
        c(0,x$dposterior(tau=x$summary["median","tau"])), col=grDevices::rgb(0.4,0.4,0.6))
  # draw density:
  lines(tau, dens, col="blue3")
  # draw prior density:
  if (prior && x$tau.prior.proper)
    lines(tau, x$dprior(tau=tau), col="red3", lty=2)
  for (i in 1:x$d) {
    xrange <- x$qposterior(beta=c(0.0005, 0.9995), which.beta=i)
    beta <- seq(xrange[1], xrange[2], le=400)
    dens <- x$dposterior(beta=beta, which.beta=i)
    plot(xrange + c(1,-1)*diff(xrange)/20, c(0, max(dens)), type="n",
         main=colnames(x$X)[i], xlab=bquote(beta[.(i)]), ylab="density")
    # draw left tail:
    idx <- which(beta < x$summary["95% lower",1+i])
    polygon(c(beta[1], beta[idx], rep(x$summary["95% lower",1+i],2)),
            c(0, dens[idx], x$dposterior(beta=x$summary["95% lower",1+i], which.beta=i), 0),
            col="grey90", border=NA)
    # draw central bulk:
    idx <- which((beta > x$summary["95% lower",1+i]) & (beta < x$summary["95% upper",1+i]))
    polygon(c(rep(x$summary["95% lower",1+i],2), beta[idx], rep(x$summary["95% upper",1+i],2)),
            c(0, x$dposterior(beta=x$summary["95% lower",1+i], which.beta=i),
              dens[idx], x$dposterior(beta=x$summary["95% upper",1+i], which.beta=i), 0),
            col="grey80", border=NA)
    # draw right tail:
    idx <- which(beta > x$summary["95% upper",1+i])
    polygon(c(rep(x$summary["95% upper",1+i],2), beta[idx], beta[length(beta)]),
            c(0, x$dposterior(beta=x$summary["95% upper",1+i], which.beta=i), dens[idx], 0),
            col="grey90", border=NA)
    # draw axes:
    abline(h=0, v=0, col="grey40")
    # draw median:
    lines(rep(x$summary["median",1+i],2),
          c(0,x$dposterior(beta=x$summary["median",1+i], which.beta=i)), col=grDevices::rgb(0.4,0.4,0.6))
    abline(h=0, v=0, col="grey40")
    # draw density:
    lines(beta, dens, col="blue3")
    # draw prior density:
    if (prior && x$beta.prior.proper[i])
      lines(beta, dnorm(x=beta, mean=x$beta.prior$mean[i],
                       sd=sqrt(x$beta.prior$covariance[i,i])),
            col="red3", lty=2)
  }
  invisible()
}



pairs.bmr <- function(x, ...)
{
  par.mfrow <- par("mfrow")
  par.mar   <- par("mar")
  par.mgp   <- par("mgp")
  on.exit(par(mfrow=par.mfrow, mar=par.mar, mgp=par.mgp))
  par(mfrow=c(x$d+1, x$d+1),
      mar=c(4,4,1,1)+0.1,
      mgp=c(2.5, 0.7, 0.0))
  diaglabels <- expression("heterogeneity ("*tau*")")
  for (i in 1:x$d)
    diaglabels <- c(diaglabels,
                    bquote(.(x$variables[i]) * " (" * beta[.(i)] * ")"))
  for (i in 1:(x$d+1)) {
    for (j in 1:(x$d+1)) {
      if (i > j) {
        plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
      } else if (i==j) {
        plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
        text(0.5, 0.5, diaglabels[i])
        box()
      } else {
        if (i==1) { # tau vs. beta[j]
          xrange <- x$qposterior(beta=c(0.0005, 0.9995), which.beta=(j-1))
          beta1  <- seq(xrange[1], xrange[2], le=50)
          yrange <- c(0, x$qposterior(tau=0.9990))
          tau    <- seq(yrange[1], yrange[2], le=50)
          xy <- expand.grid(beta1, tau)
          # evaluate posterior density at grid points:
          post <- matrix(x$dposterior(tau=xy[,2], beta=xy[,1], which.beta=(j-1), log=TRUE),
                         nrow=length(beta1), ncol=length(tau))
          ## determine MAP value:
          map.value <- max(post[is.finite(post)])
          image(beta1, tau, exp(post), axes=FALSE,
                col=grey((seq(1,0,le=128))^2),
                breaks=seq(0,exp(map.value),length=129),
                xlab=bquote(beta[.(j-1)]), ylab=expression(tau), main="")
          abline(h=0, v=0, col="grey40")
          contour(beta1, tau, post-map.value, add=TRUE, col="red",
                  levels=-0.5*qchisq(p=c(0.5, 0.9, 0.95, 0.99), df=2),
                  labels=paste(c(50, 90, 95, 99),"%",sep=""))
          
          axis(1); axis(2); box()
        } else {    # beta[i] vs. beta[j]
          xrange <- x$qposterior(beta=c(0.0005, 0.9995), which.beta=(j-1))
          beta1  <- seq(xrange[1], xrange[2], le=50)
          yrange <- x$qposterior(beta=c(0.0005, 0.9995), which.beta=(i-1))
          beta2  <- seq(yrange[1], yrange[2], le=50)
          xy <- expand.grid(beta1, beta2)
          # evaluate posterior density at grid points:
          post <- matrix(x$dposterior(beta=xy, which.beta=c(j-1, i-1), log=TRUE),
                         nrow=length(beta1), ncol=length(beta2))
          ## determine (bivariate) MAP value:
          map.value <- max(post[is.finite(post)])
          opti <- optim(x$MAP[2,c(j,i)],
                        function(xx){-x$dposterior(beta=xx, which.beta=c(j-1,i-1), log=TRUE)})
          map.arg   <- opti$par
          map.value <- -opti$value

          image(beta1, beta2, exp(post), axes=FALSE,
                col=grey((seq(1,0,le=128))^2),
                breaks=seq(0,exp(map.value),length=129),
                xlab=bquote(beta[.(j-1)]), ylab=bquote(beta[.(i-1)]), main="")
          abline(h=0, v=0, col="grey40")
          contour(beta1, beta2, post-map.value, add=TRUE, col="red",
                  levels=-0.5*qchisq(p=c(0.5, 0.9, 0.95, 0.99), df=2),
                  labels=paste(c(50, 90, 95, 99),"%",sep=""))
          points(map.arg[1], map.arg[2], col="red", pch=3)

          axis(1); axis(2); box()
        }
      }
    }
  }
  invisible()
}



forestplot.bmr <- function(x,
                           X.mean, X.prediction,
                           labeltext,
                           exponentiate  = FALSE,
                           shrinkage     = TRUE,
                           heterogeneity = TRUE,
                           digits        = 2,
                           decplaces.X,
                           plot          = TRUE,
                           fn.ci_norm, fn.ci_sum, col, legend=NULL, boxsize, ...)
# ARGUMENTS:
#   x            :  a "bmr" object.
#   X.mean       :  regressor matrix for summary (mean) estimate(s)
#   X.prediction :  regressor metrix for predictions
#   labeltext    :  you may provide an alternative "labeltext" argument here
#                   (see also the "forestplot()" help).
#   exponentiate :  flag indicating whether to exponentiate numbers (figure and table).
#   shrinkage    :  flag indicating whether to show shrinkage estimates.
#   digits       :  number of significant digits to be shown (based on standard deviations).
#   decplaces.X  :  number of decimal places to show for regressor matrices 
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
  # auxiliary functions:
  count.decimals <- function(x, max=3)
  # try to count relevant number of decimal places (after decimal point)
  {
    stopifnot(is.vector(x), is.numeric(x), max >= 0)
    i <- 0
    while ((i<max) && (!all((x-round(x,i)) < sqrt(.Machine$double.eps)))) {
      i <- i+1
    }
    return(i)
  }
  decplaces <- function(x, signifdigits=3)
  # number of decimal places (after decimal point)
  # to be displayed if you want at least "signifdigits" digits
  {
    return(max(c(0, -(floor(log10(x))-(signifdigits-1)))))
  }
  # some sanity checks for the provided arguments:
  stopifnot(is.element("bmr", class(x)),
            length(digits)==1, digits==round(digits), digits>=0,
            length(exponentiate)==1, is.logical(exponentiate),
            length(shrinkage)==1, is.logical(shrinkage),
            length(plot)==1, is.logical(plot))
  #########################################
  # default treatment of "X.mean" argument:
  if (missing(X.mean)) {
    X.mean <- diag(1, nrow=x$d, ncol=x$d)
    dimnames(X.mean) <- list(x$variables, x$variables)
    meanNumber <- x$d
  } else {
    if (is.null(X.mean) || all(is.na(X.mean))) {
      X.mean <- NULL
      meanNumber <- 0
    } else {
      stopifnot(is.numeric(X.mean),
                (is.matrix(X.mean) && (ncol(X.mean) == x$d))
                | is.vector(X.mean) && (length(X.mean) == x$d))
      if (is.vector(X.mean)) {
        X.mean <- matrix(X.mean, nrow=1, ncol=x$d)
      }
      if (is.null(colnames(X.mean))) colnames(X.mean) <- x$variables
      meanNumber <- nrow(X.mean)
      if (is.null(rownames(X.mean))) rownames(X.mean) <- sprintf("mean%02d", 1:meanNumber)
    }
  }
  # default treatment of "X.prediction" argument:
  if (missing(X.prediction) || (is.null(X.prediction) || all(is.na(X.prediction)))) {
      X.prediction <- NULL
      predNumber <- 0
  } else {
    stopifnot(is.numeric(X.prediction),
              (is.matrix(X.prediction) && (ncol(X.prediction) == x$d))
              | (is.vector(X.prediction) && (length(X.prediction) == x$d)))
    if (is.vector(X.prediction)) {
      X.prediction <- matrix(X.prediction, nrow=1, ncol=x$d)
    }
    if (is.null(colnames(X.prediction))) colnames(X.prediction) <- x$variables
    predNumber <- nrow(X.prediction)
    if (is.null(rownames(X.prediction))) rownames(X.prediction) <- sprintf("prediction%02d", 1:predNumber)
  }
  ############################################
  # plotting data (1) -- the quoted estimates:
  q95 <- qnorm(0.975)
  ma.dat <- rbind(NA_real_,
                  cbind(x$y, x$y - q95*x$sigma, x$y + q95*x$sigma),
                  matrix(NA_real_, nrow=(meanNumber + predNumber), ncol=3))
  colnames(ma.dat) <- c("estimate", "lower", "upper")
  rownames(ma.dat) <- c("", x$label, rep("", meanNumber + predNumber))
  if (meanNumber > 0)  { # compute summaries:
    rownames(ma.dat)[1 + x$k + (1:meanNumber)] <- rownames(X.mean)
    for (i in 1:meanNumber) {
      ma.dat[1+x$k+i, "estimate"] <- x$qpredict(p=0.5, x=X.mean[i,], mean=TRUE)
      ma.dat[1+x$k+i, c("lower","upper")] <- x$pred.interval(x=X.mean[i,],
                                                             mean=TRUE,
                                                             method=x$interval.type)
    }
  }
  if (predNumber > 0) { # compute predictions:
    rownames(ma.dat)[1 + x$k + meanNumber + (1:predNumber)] <- rownames(X.prediction)
    for (i in 1:predNumber) {
      ma.dat[1+x$k+meanNumber+i, "estimate"] <- x$qpredict(p=0.5, x=X.prediction[i,], mean=FALSE)
      ma.dat[1+x$k+meanNumber+i, c("lower","upper")] <- x$pred.interval(x=X.prediction[i,],
                                                                       mean=FALSE,
                                                                       method=x$interval.type)
    }
  }
  # plotting data (2) -- the shrinkage estimates:
  ma.shrink <- rbind(NA_real_,
                     t(x$theta)[,c("median","95% lower","95% upper")],
                     matrix(NA_real_, nrow=(meanNumber + predNumber), ncol=3))
  colnames(ma.shrink) <- c("estimate", "lower", "upper")
  rownames(ma.shrink) <- c("", x$label, rep("", meanNumber + predNumber))
  if (meanNumber > 0)  rownames(ma.shrink)[1 + x$k + (1:meanNumber)] <- rownames(X.mean)
  if (predNumber > 0) rownames(ma.shrink)[1 + x$k + meanNumber + (1:predNumber)] <- rownames(X.prediction)
  if (exponentiate) {
    ma.dat    <- exp(ma.dat)
    ma.shrink <- exp(ma.shrink)
  }
  # generate "labeltext" data table for plot (unless already provided):
  if (missing(labeltext)) {
    # determine numbers of digits based on standard deviations:
    stdevs <- apply(ma.dat[,c("lower", "upper")], 1, diff) / (2*q95) # (simple proxy for std.devs.) 
    stdevs <- abs(stdevs[is.finite(stdevs) & (stdevs != 0)])
    formatstring <- paste0("%.", decplaces(stdevs, digits), "f")
    # set up data table:
    labeltext <- matrix(NA_character_, nrow=nrow(ma.dat), ncol=3+x$d)
    # 1st row & column:
    labeltext[1,] <- c("study", x$variables, "estimate", "95% CI")
    labeltext[1+(1:x$k),1] <- x$labels
    if (meanNumber > 0) labeltext[1 + x$k + (1:meanNumber), 1] <- rownames(X.mean)
    if (predNumber > 0) labeltext[1 + x$k + meanNumber + (1:predNumber), 1] <- rownames(X.prediction)
    # regressor matrix; determine number of decimal places to show:
    if (missing(decplaces.X)) {
      decplaces.X <- rep(NA_integer_, x$d)
      for (i in 1:x$d) decplaces.X[i] <- count.decimals(rbind(x$X, X.mean, X.prediction)[,i])
    } else {
      stopifnot(is.vector(decplaces.X), all(is.finite(decplaces.X)),
                all(decplaces.X >= 0), all(decplaces.X==round(decplaces.X)))
      if (length(decplaces.X) != x$d) decplaces.X <- rep(decplaces.X[1], x$d)
    }
    names(decplaces.X) <- x$variables
    # regressor matrix; convert to text:
    for (i in 1:x$d) {
      labeltext[1+(1:x$k),1+i] <- sprintf(paste0("%.", decplaces.X[i], "f"), x$X[,i])
      if (meanNumber > 0)
        labeltext[1+x$k+(1:meanNumber),1+i] <- sprintf(paste0("%.", decplaces.X[i], "f"), X.mean[,i])
      if (predNumber > 0)
        labeltext[1+x$k+meanNumber+(1:predNumber),1+i] <- sprintf(paste0("%.", decplaces.X[i], "f"), X.prediction[,i])
    }
    # estimates & CIs:
    for (i in 2:(nrow(ma.dat))) {
      labeltext[i,x$d+2] <- sprintf(formatstring, ma.dat[i,"estimate"])
      labeltext[i,x$d+3] <- paste0("[", sprintf(formatstring, ma.dat[i,"lower"]),
                                   ", ", sprintf(formatstring, ma.dat[i,"upper"]), "]")
    }
  }
  # add horizontal lines to plot:
  if ((meanNumber>0) & (predNumber>0)) {         # (both summaries & predictions)
    horizl <- list(grid::gpar(col="grey"), grid::gpar(col="grey"), grid::gpar(col="grey"))
    names(horizl) <- as.character(c(2, 2+x$k, 2+x$k+meanNumber))
  } else if ((meanNumber>0) | (predNumber>0)) {  # (only one but not both)
    horizl <- list(grid::gpar(col="grey"), grid::gpar(col="grey"))
    names(horizl) <- as.character(c(2, 2+x$k))
  } else {                                      # (studies only)
    horizl <- list(grid::gpar(col="grey"))
    names(horizl) <- as.character(2)
  }
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
    for (i in 1:(1+x$k+meanNumber))
      fn.ci_sum[[i]] <- function(y.offset,...) {forestplot::fpDrawSummaryCI(y.offset=0.5,...)}
    if (predNumber>0)
      for (i in (1+x$k+meanNumber + (1:predNumber)))
        fn.ci_sum[[i]] <- function(y.offset,...) {forestplot::fpDrawBarCI(y.offset=0.5,...)}
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
    boxsize <- c(rep(0.25,x$k+1), rep(0.4, meanNumber), rep(0.2, predNumber))
  }
  if (shrinkage && is.null(legend))
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
                                 is.summary = c(TRUE, rep(FALSE, x$k), rep(TRUE,meanNumber+predNumber)),
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
  invisible(list("data"         = ma.dat[-1,],
                 "X.mean"       = X.mean,
                 "X.prediction" = X.prediction,
                 "shrinkage"    = ma.shrink[1+(1:x$k),],
                 "labeltext"    = labeltext,
                 "forestplot"   = fp))
}



traceplot.bmr <- function(x, mulim, taulim, ci=FALSE,
                          ylab="effect",
                          prior=FALSE,
                          infinity=FALSE,
                          rightmargin=8,
                          col=rainbow(x$k), labcol=col,
                          X, Xlabels,
                          Xcols="black", Xlabcols=Xcols,
                          ...)
{
  stopifnot(missing(mulim) || (length(mulim) == 2),
            missing(taulim) || (length(taulim) <= 2),
            is.character(ylab), length(ylab)==1,
            is.logical(prior), length(prior)==1,
            is.logical(infinity), length(infinity)==1,
            rightmargin >= 0, ((length(col)==x$k) | (length(col)==1)))
  q975 <- qnorm(0.975)
  # specify colors for axes etc.:
  colvec <- c("axis"   = "grey40",
              "grid"   = "grey85",
              "median" = "grey60",
              "ci"     = "grey80",
              "tail"   = "grey90")
  if (length(col)==1) col <- rep(col, x$k)
  if (infinity & any(x$beta.prior.proper)) {
    warning("beta prior ignored for `tau=Inf' computations!")
  }
  if (prior & (!x$tau.prior.proper)) {
    warning("Note that plots of improper priors may not be sensibly scaled.")
  }
  # default treatment of "X" argument:
  if (missing(X) || (is.null(X) || all(is.na(X)))) {
      X <- NULL
      lincoNumber <- 0
  } else {
    stopifnot(is.numeric(X),
              (is.matrix(X) && (ncol(X) == x$d))
              || (is.vector(X) && (length(X) == x$d)))
    if (is.vector(X)) {
      X <- matrix(X, nrow=1, ncol=x$d)
    }
    if (is.null(colnames(X))) colnames(X) <- x$variables      
    lincoNumber <- nrow(X)
  }
  # determine "X" labels and colours:
  if (lincoNumber > 0) {
    if (!missing(Xlabels)) {
      stopifnot(is.character(Xlabels), length(Xlabels) == nrow(X))
      rownames(X) <- Xlabels
    } else {
      Xlabels <- rownames(X)   # take from X row names...
      if (is.null(Xlabels)) {  # ...or make up some
        Xlabels <- sprintf("X.%02d", 1:nrow(X))
      }
    }
    stopifnot(is.vector(Xcols), ((length(Xcols)==1) | (length(Xcols)==lincoNumber)))
    if (length(Xcols)==1) Xcols <- rep(Xcols, lincoNumber)
  }
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
  
  if (infinity) {
    xlim <- taurange + c(0, 0.15) * diff(taurange)
    infx <- xlim[2] + 0.04*diff(xlim) # the "infinity" x-coordinate
  } else {
    xlim <- taurange
    infx <- NA_real_
  }
  
  vertlines <- pretty(taurange)
  # ensure no tickmarks beyond plotted tau range:
  if (max(vertlines) > (taurange[2] + 0.04*diff(taurange)))
    vertlines <- vertlines[-length(vertlines)]
  
  if (missing(mulim)) mulim <- NULL
  
  mutrace <- function(x)
  {
    # vector of tau values:
    tau <- seq(max(c(0,taurange[1]-0.04*diff(taurange))),
                   taurange[2]+0.04*diff(taurange), le=200)
    # moments for individual studies:
    cm.indiv <- array(NA_real_, dim=c(length(tau), 2, x$k),
                      dimnames=list(NULL, c("mean", "sd"), x$labels))
    for (i in 1:x$k) {
      cm.indiv[,,i] <- x$shrink.moment(tau=tau, which=i)
    }
    # moments for linear combinations (if applicable):
    if (lincoNumber > 0) {
      cm.linco <- array(NA_real_, dim=c(length(tau), 2, lincoNumber),
                        dimnames=list(NULL, c("mean", "sd"), Xlabels))
      for (i in 1:lincoNumber) {
        cm.linco[,,i] <- x$pred.moments(tau=tau, x=X[i,])
      }
    }

    # determine axis range for "effect" (y-) axis
    if (!is.null(mulim) && (all(is.finite(mulim)) && (mulim[1] < mulim[2]))) {
      # user-defined:
      murange <- mulim
    } else {
      # based on data:
      if (ci){
        murange <- range(c(range(cm.indiv[,"mean",]-q975*cm.indiv[,"sd",]),
                           range(cm.indiv[,"mean",]+q975*cm.indiv[,"sd",])))
        if (lincoNumber > 0) {
          murange <- range(c(murange,
                             range(cm.linco[,"mean",]-q975*cm.linco[,"sd",]),
                             range(cm.linco[,"mean",]+q975*cm.linco[,"sd",])))
        }
      } else {
        murange <- range(cm.indiv[,"mean",])
        if (lincoNumber > 0) {
          murange <- range(c(murange,
                             range(cm.linco[,"mean",])))
        }
      }
      # ensure that estimates are included:
      if (infinity) murange <- range(murange, x$y)
    }
  
    plot(taurange, murange, xlim=xlim,
         type="n", axes=FALSE, xlab="", ylab=ylab, main="", ...)
    abline(v=vertlines, col=colvec["grid"])
    abline(h=pretty(murange), col=colvec["grid"])
    abline(v=0, col=colvec["axis"])
    # grey CI shading:
    if (ci) {
      for (i in 1:x$k) {
        polygon(c(tau, rev(tau)),
                c(cm.indiv[,"mean",i] - q975*cm.indiv[,"sd",i],
                  rev(cm.indiv[,"mean",i] + q975*cm.indiv[,"sd",i])),
                col=grey(0.75, alpha=0.25), border=NA)
      }
      if (lincoNumber > 0) {
        for (i in 1:lincoNumber) {
          polygon(c(tau, rev(tau)),
                  c(cm.linco[,"mean",i] - q975*cm.linco[,"sd",i],
                    rev(cm.linco[,"mean",i] + q975*cm.linco[,"sd",i])),
                  col=grey(0.75, alpha=0.25), border=NA)
        }
      }
    }    
    # individual estimates:
    matlines(tau, cm.indiv[,"mean",], col=col, lty=1)
    if (ci) {
      matlines(tau, cm.indiv[,"mean",]-q975*cm.indiv[,"sd",], col=col, lty=3)
      matlines(tau, cm.indiv[,"mean",]+q975*cm.indiv[,"sd",], col=col, lty=3)
    }
    # linear combinations:
    if (lincoNumber > 0) {
      matlines(tau, cm.linco[,"mean",], col=Xcols, lty=2, lwd=1.5)
      if (ci) {
        matlines(tau, cm.linco[,"mean",]-q975*cm.linco[,"sd",],
                 col=Xcols, lty=3, lwd=1.5)
        matlines(tau, cm.linco[,"mean",]+q975*cm.linco[,"sd",],
                 col=Xcols, lty=3, lwd=1.5)
      }
    }
      
    if (infinity) {
      # individual studies:
      labpos.indiv   <- x$y
      for (i in 1:x$k) {
        lines(c(max(tau), infx),
              c(cm.indiv[length(tau),"mean",i], labpos.indiv[i]),
              col=col[i], lty="13", lwd=1.5)
      }
      # linear combinations:
      if (lincoNumber > 0) {
        # compute (unweighted) least-squares fit:
        lscoef <- stats::lm.fit(x=x$X, y=x$y)$coefficients
        labpos.linco <- X %*% lscoef
        for (i in 1:lincoNumber) {
          lines(c(max(tau), infx),
                c(cm.linco[length(tau),"mean",i], labpos.linco[i]),
                col=Xcols[i], lty="13", lwd=2)
        }
      }
    } else {
      labpos.indiv <- cm.indiv[length(tau),"mean",]
      if (lincoNumber > 0) {
        labpos.linco <- cm.linco[length(tau),"mean",]
      }
    }
    axis(2)
    for (i in 1:x$k)
      axis(side=4, at=labpos.indiv[i],
           labels=x$labels[i], tick=FALSE,
           col.axis=labcol[i], las=1)
    if (lincoNumber > 0) {
      for (i in 1:lincoNumber) {
        axis(side=4, at=labpos.linco[i],
             labels=Xlabels[i], tick=FALSE,
             col.axis=Xlabcols[i], las=1)
      }
    }
    invisible()
  }
  
  taumarginal <- function(x)
  # NB: function is (essentially) identical to the one within "plot.bayesmeta()"
  {
    # range of tau values:
    tau <- seq(max(c(0,taurange[1]-0.04*diff(taurange))),
                   taurange[2]+0.04*diff(taurange), le=200)
    # corresponding posterior density:
    dens <- x$dposterior(tau=tau)
    # empty plot:
    maxdens <- max(dens[is.finite(dens)],na.rm=TRUE)
    plot(c(taurange[1],taurange[2]), c(0,maxdens), xlim=xlim,
         type="n", axes=FALSE, xlab="", ylab="", main="")
    abline(v=vertlines, col=colvec["grid"])
    # "fix" diverging density:
    dens[!is.finite(dens)] <- 10*maxdens
    # light grey shaded contour for density across whole range:
    polygon(c(0,tau,max(tau)), c(0,dens,0), border=NA, col=colvec["tail"])
    # dark grey shaded contour for density within 95% bounds:
    indi <- ((tau>=x$summary["95% lower","tau"]) & (tau<=x$summary["95% upper","tau"]))
    polygon(c(rep(x$summary["95% lower","tau"],2), tau[indi], rep(x$summary["95% upper","tau"],2)),
            c(0, min(c(x$dposterior(tau=x$summary["95% lower","tau"]), 10*maxdens)),
              dens[indi], x$dposterior(tau=x$summary["95% upper","tau"]), 0),
            border=NA, col=colvec["ci"])
    # vertical line at posterior median:
    lines(rep(x$summary["median","tau"],2), c(0,x$dposterior(tau=x$summary["median","tau"])),
          col=colvec["median"])
    # actual density line:
    lines(tau, dens, col="black")
    # y-axis:
    abline(v=0, col=colvec["axis"])
    # x-axis:
    lines(taurange + c(-1,1) * 0.04*diff(taurange), c(0,0), col=colvec["axis"])
    # plot prior density (if requested):
    if (prior) {
      lines(tau, x$dprior(tau=tau), col="black", lty="dashed")
    }
    # add axes, labels, bounding box, ...
    mtext(side=1, line=par("mgp")[1], expression("heterogeneity "*tau))
    #mtext(side=2, line=par("mgp")[2], expression("marginal posterior density"))
    if (infinity) {
      axis(1, at=c(vertlines, infx),
           labels=c(as.numeric(vertlines), expression(infinity)))
    } else {
      axis(1, at=vertlines)
    }
    invisible()
  }
  
  # make sure to properly re-set graphical parameters later:
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
  invisible(X)
}
