\name{Cochran1954}
\docType{data}
\alias{Cochran1954}
\title{Fly counts example data}
\description{This data set gives average estimated counts of flies along
  with standard errors from 7 different observers.}
\usage{data("Cochran1954")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{observer} \tab \code{character} \tab identifier \cr
    \bold{mean}     \tab \code{numeric}   \tab mean count \cr
    \bold{se2}      \tab \code{numeric}   \tab \emph{squared} standard error 
  }
}
\details{
  Quoting from Cochran (1954), example 3, p.119:
  \dQuote{In studies by the U.S. Public Health Service of observers'
    abilities to count the number of flies which settle momentarily on a
    grill, each of 7 observers was shown, for a brief period, grills
    with known numbers of flies impaled on them and asked to estimate
    the numbers. For a given grill, each observer made 5 independent
    estimates. The data in table 9 are for a grill which actually
    contained 161 flies. Estimated variances are based on 4 degrees of
    freedom each. [...] The only point of interest in estimating the
    overall mean is to test whether there is any consistent bias among
    observers in estimating the 161 flies on the grill. Although
    inspection of table 9 suggests no such bias, the data will serve
    to illustrate the application of partial weighting.} 
}
\source{W.G. Cochran. \href{http://dx.doi.org/10.2307/3001666}{The
    combination of estimates from different experiments}.
    \emph{Biometrics}, \bold{10}(1):101-129, 1954.}
\examples{
\dontrun{
data("Cochran1954")
# analysis using improper uniform prior
# (may take a few seconds to compute!):
bma <- bayesmeta(y=Cochran1954[,"mean"], sigma=sqrt(Cochran1954[,"se2"]),
                 label=Cochran1954[,"observer"])

# show joint posterior density:
plot(bma, which=2, main="Cochran example")
# show (known) true parameter value:
abline(h=161)
}
}
\keyword{datasets}