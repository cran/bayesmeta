\name{SidikJonkman2007}
\docType{data}
\alias{SidikJonkman2007}
\title{Postoperative complication odds example data}
\description{This data set contains the outcomes from 29 randomized
  clinical trials comparing the odds of postoperative complications in
  laparoscopic inguinal hernia repair (LIHR) versus conventional open
  inguinal hernia repair (OIHR).}
\usage{data("SidikJonkman2007")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{id}          \tab \code{character} \tab identifier used in 
    original publication by Memon et al. (2003)\cr
    \bold{id.sj}       \tab \code{numeric}   \tab identifier used by 
    Sidik and Jonkman (2007)\cr
    \bold{year}        \tab \code{numeric}   \tab publication year\cr
    \bold{lihr.events} \tab \code{numeric}   \tab number of events under LIHR\cr
    \bold{lihr.cases}  \tab \code{numeric}   \tab number of cases under LIHR\cr
    \bold{oihr.events} \tab \code{numeric}   \tab number of events under OIHR\cr
    \bold{oihr.cases}  \tab \code{numeric}   \tab number of cases under OIHR
  }
}
\details{Analysis may be done based on the logarithmic odds ratios:

  log(\code{lihr.events}) - log(\code{lihr.cases}-\code{lihr.events}) -
  log(\code{oihr.events}) + log(\code{oihr.cases}-\code{oihr.events})

  and corresponding standard errors:

  sqrt(1/\code{lihr.events} + 1/(\code{lihr.cases}-\code{lihr.events}))
       + 1/\code{oihr.events} + 1/(\code{oihr.cases}-\code{oihr.events}))

  (you may also leave these computations to the \pkg{metafor} package's
  \code{\link[metafor]{escalc}()} function).
       
  The data set was used to compare different estimators for the
  (squared) heterogeneity \eqn{\tau^2}. The values yielded for this data
  set were (see Tab.1 in Sidik and Jonkman (2007)):
  \tabular{ll}{
    method of moments (MM)    \tab 0.429 \cr
    variance component (VC)   \tab 0.841 \cr
    maximum likelihood (ML)   \tab 0.562 \cr
    restricted ML (REML)      \tab 0.598 \cr
    empirical Bayes (EB)      \tab 0.703 \cr
    model error variance (MV) \tab 0.818 \cr
    variation of MV (MVvc)    \tab 0.747 
  }
}
\source{M.A. Memon, N.J. Cooper, B. Memon, M.I. Memon, and
  K.R. Abrams.
  Meta-analysis of randomized clinical trials comparing open and
  laparoscopic inguinal hernia repair. 
  \emph{British Journal of Surgery}, \bold{90}(12):1479-1492, 2003.
  \doi{10.1002/bjs.4301}.
}
\references{K. Sidik and J.N. Jonkman.
  A comparison of heterogeneity variance estimators in combining results
  of studies.
  \emph{Statistics in Medicine}, \bold{26}(9):1964-1981, 2007.
  \doi{10.1002/sim.2688}.
}
\examples{
data("SidikJonkman2007")
# add log-odds-ratios and corresponding standard errors:
sj <- SidikJonkman2007
sj <- cbind(sj, "log.or"=log(sj[,"lihr.events"])-log(sj[,"lihr.cases"]-sj[,"lihr.events"])
                             -log(sj[,"oihr.events"])+log(sj[,"oihr.cases"]-sj[,"oihr.events"]),
            "log.or.se"=sqrt(1/sj[,"lihr.events"] + 1/(sj[,"lihr.cases"]-sj[,"lihr.events"])
                             + 1/sj[,"oihr.events"] + 1/(sj[,"oihr.cases"]-sj[,"oihr.events"])))

\dontrun{
# analysis using weakly informative Cauchy prior
# (may take a few seconds to compute!):
ma <- bayesmeta(y=sj[,"log.or"], sigma=sj[,"log.or.se"], label=sj[,"id.sj"],
                tau.prior=function(t){dhalfcauchy(t,scale=1)})

# show heterogeneity's posterior density:
plot(ma, which=4, main="Sidik/Jonkman example", prior=TRUE)

# show some numbers (mode, median and mean):
abline(v=ma$summary[c("mode","median","mean"),"tau"], col="blue")

# compare with Sidik and Jonkman's estimates:
sj.estimates <- sqrt(c("MM"  = 0.429,   # method of moments estimator
                       "VC"  = 0.841,   # variance component type estimator
                       "ML"  = 0.562,   # maximum likelihood estimator
                       "REML"= 0.598,   # restricted maximum likelihood estimator
                       "EB"  = 0.703,   # empirical Bayes estimator
                       "MV"  = 0.818,   # model error variance estimator
                       "MVvc"= 0.747))  # a variation of the MV estimator
abline(v=sj.estimates, col="red", lty="dashed")

# generate forest plot:
fp <- forestplot(ma, exponentiate=TRUE, plot=FALSE)
# add extra columns for ID and year:
labtext <- fp$labeltext
labtext[1,1] <- "ID 2"
labtext[31:32,1] <- ""
labtext <- cbind(c("ID 1", SidikJonkman2007[,"id"], "mean","prediction"),
                 labtext[,1],
                 c("year", as.character(SidikJonkman2007[,"year"]), "", ""),
                 labtext[,-1])
# plot:
forestplot(ma, labeltext=labtext, exponentiate=TRUE, 
           xlog=TRUE, xlab="odds ratio", xticks=c(0.1,1,10))

}
}
\keyword{datasets}