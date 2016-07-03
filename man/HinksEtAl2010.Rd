\name{HinksEtAl2010}
\docType{data}
\alias{HinksEtAl2010}
\title{JIA example data}
\description{Log odds ratios indicating association of a genetic variant
  (CCR5) with juvenile idiopathic arthritis (JIA).}
\usage{data("HinksEtAl2010")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{study}     \tab \code{character} \tab publication identifier \cr
    \bold{year}      \tab \code{numeric}   \tab publication year \cr
    \bold{country}   \tab \code{character} \tab country \cr
    \bold{or}        \tab \code{numeric}   \tab odds ratio (OR)\cr
    \bold{or.lower}  \tab \code{numeric}   \tab lower 95 percent confidence bound for OR \cr
    \bold{or.upper}  \tab \code{numeric}   \tab upper 95 percent confidence bound for OR \cr
    \bold{log.or}    \tab \code{numeric}   \tab logarithmic OR \cr
    \bold{log.or.se} \tab \code{numeric}   \tab standard error of logarithmic OR \cr
  }
}
\details{
  Results from a genetic association study (Hinks et al, 2010) were
  combined with data from two additional studies (Prahalad et al., 2006;
  Lindner et al., 2007) in order to determine the combined evidence
  regarding the association of a particular genetic marker (CCR5) with
  juvenile idiopathic arthritis (JIA).
}
\source{
  A. Hinks et al.
  \href{http://dx.doi.org/10.1038/gene.2010.25}{Association of the CCR5
    gene with juvenile idiopathic arthritis}. \emph{Genes and Immunity},
  \bold{11}(7):584-589, 2010. 
}
\references{
  S. Prahalad et al.
  \href{http://dx.doi.org/10.1038/sj.gene.6364317}{Association of two
    functional polymorphisms in the CCR5 gene with juvenile rheumatoid
    arthritis}
  \emph{Genes and Immunity}, \bold{7}:468-475, 2006.

  E. Lindner et al.
  \href{http://dx.doi.org/10.1186/1471-2350-8-33}{Lack of association
  between the chemokine receptor 5 polymorphism CCR5delta32 in
  rheumatoid arthritis and juvenile idiopathic arthritis}.
  \emph{BMC Medical Genetics}, \bold{8}:33, 2007.

  C. Roever, G. Knapp, T. Friede.
  \href{http://dx.doi.org/10.1186/s12874-015-0091-1}{Hartung-Knapp-Sidik-Jonkman
  approach and its modification for random-effects meta-analysis with
  few studies}. \emph{BMC Medical Research Methodology}, 15:99, 2015.
}
\examples{
\dontrun{
data("HinksEtAl2010")

# perform meta analysis based on weakly informative half-normal prior:
bma01 <- bayesmeta(y=HinksEtAl2010$log.or,
                   sigma=HinksEtAl2010$log.or.se,
                   lab=HinksEtAl2010$study,
                   tau.prior=function(t){dhalfnormal(t,scale=1.0)})

# perform meta analysis based on slightly more informative half-normal prior:
bma02 <- bayesmeta(y=HinksEtAl2010$log.or,
                   sigma=HinksEtAl2010$log.or.se,
                   lab=HinksEtAl2010$study,
                   tau.prior=function(t){dhalfnormal(t,scale=0.5)})

# show heterogeneity posteriors:
plot(bma01, which=4, prior=TRUE)
plot(bma02, which=4, prior=TRUE)

# show heterogeneity estimates:
rbind("half-normal(1.0)"=bma01$summary[,"tau"],
      "half-normal(0.5)"=bma02$summary[,"tau"])
# show q-profile confidence interval for tau in comparison:
require("metafor")
ma03 <- rma.uni(yi=log.or, sei=log.or.se, slab=study, data=HinksEtAl2010)
confint(ma03)$random["tau",c("ci.lb","ci.ub")]
# show I2 values in the relevant range:
tau <- seq(0, 0.7, by=0.1)
cbind("tau"=tau,
      "I2" =bma01$I2(tau=tau))

# show effect estimates:
round(rbind("half-normal(1.0)" = bma01$summary[,"mu"],
            "half-normal(0.5)" = bma02$summary[,"mu"]), 5)

# show forest plot:
forest(bma02)
# show the "shrunk" estimate of the Hinks (2006) study in comparison:
points(bma02$theta["median","Hinks (2006)"], -3.2, pch=22, col="red")
lines(bma02$theta[c("95% lower","95% upper"),"Hinks (2006)"], rep(-3.2,2), col="red")
}
}
\keyword{datasets}