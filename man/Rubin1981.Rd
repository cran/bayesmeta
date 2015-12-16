\name{Rubin1981}
\docType{data}
\alias{Rubin1981}
\title{8-schools example data}
\description{SAT coaching experiments in 8 schools.}
\usage{data("Rubin1981")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{school} \tab \code{character} \tab school identifier \cr
    \bold{effect} \tab \code{numeric}   \tab effect estimate \cr
    \bold{stderr} \tab \code{numeric}   \tab associated standard error 
  }
}
\details{
  Quoting from Gelman et al. (1997), Sec. 5.5:
  \dQuote{A study was performed for the Educational Testing Service to
  analyze the effects of special coaching programs for SAT-V
  (Scholastic Aptitude Test-Verbal) in each of eight high schools. The
  outcome variable in each study was the score on a special
  administration of the SAT-V, a standardized multiple choice test
  administered by the Educational Testing Service and used to help
  colleges make admissions decisions; the scores can vary between 200
  and 800, with mean about 500 and standard deviation about 100. The SAT
  examinations are designed to be resistant to short-term efforts
  directed specifically toward improving performance on the test;
  instead they are designed to reflect knowledge acquired and abilities
  developed over many years of education. Nevertheless, each of the
  eight schools in this study considered its short-term coaching program
  to be very successful at increasing SAT scores. Also, there was no
  prior reason to believe that any of the eight programs was more
  effective than any other or that some were more similar in effect to
  each other than to any other.}
}
\source{A. Gelman, J.B. Carlin, H. Stern, and D.B. Rubin.
  \emph{Bayesian data analysis}. Chapman & Hall / CRC, Boca Raton, 1997.} 
\references{
  D.B. Rubin. \href{http://dx.doi.org/10.3102/10769986006004377}{Estimation
    in parallel randomized experiments}.
  \emph{Journal of Educational Statistics}, \bold{6}(4):377-401, 1981.

  A. Gelman. Prior distributions for variance parameters in hierarchical
  models. \emph{Bayesian Analysis}, \bold{1}(3):515-534, 2006.
}
\examples{
\dontrun{
data(Rubin1981)
# analysis using weakly informative half-Cauchy prior
# (may take a few seconds to compute!):
schools <- bayesmeta(y=Rubin1981[,"effect"], sigma=Rubin1981[,"stderr"],
                     labels=Rubin1981[,"school"],
                     tau.prior=function(x){return(dhalfcauchy(x, scale=25))})

# show summary:
summary(schools)
}
}
\keyword{datasets}