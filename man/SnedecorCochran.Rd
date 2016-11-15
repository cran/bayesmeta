\name{SnedecorCochran}
\docType{data}
\alias{SnedecorCochran}
\title{Artificial insemination of cows example data}
\description{This data set gives means and (squared) standard errors of
  percentages of conceptions obtained from samples for six bulls.}
\usage{data("SnedecorCochran")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{no}   \tab \code{character} \tab identifier \cr
    \bold{n}    \tab \code{numeric}   \tab sample size \cr
    \bold{mean} \tab \code{numeric}   \tab mean \cr
    \bold{var}  \tab \code{numeric}   \tab variance (squared standard error)
  }
}
\details{
  Quoting from Snedecor and Cochran (1967), Sec. 10.18:
  \dQuote{In research on artificial insemination of cows, a series of
    semen samples from a bull are sent out and tested for their ability
    to produce conceptions. The following data from a larger set kindly
    supplied by Dr. G. W. Salisbury, show the percentages of conceptions
    obtained from the samples for six bulls.}
}
\source{J. Hartung, G. Knapp, and B.K. Sinha. \emph{Statistical
  meta-analysis with applications}. Wiley, Hoboken, NJ, USA, 2008.}
\references{G.W. Snedecor and W.G. Cochran. \emph{Statistical Methods}. 
  Iowa State University Press, Ames, IA, USA, 6th edition, 1967.}
\examples{
data("SnedecorCochran")

\dontrun{
# analyze using uniform prior:
bma1 <- bayesmeta(y=SnedecorCochran[,"mean"],
                  sigma=sqrt(SnedecorCochran[,"var"]),
                  label=SnedecorCochran[,"no"],
                  tau.prior="uniform")

# analyze using Jeffreys prior:
bma2 <- bayesmeta(y=SnedecorCochran[,"mean"],
                  sigma=sqrt(SnedecorCochran[,"var"]),
                  label=SnedecorCochran[,"no"],
                  tau.prior="Jeffreys")

# compare results:
print(bma1)
print(bma2)

forestplot(bma1)
forestplot(bma2)
}
}
\keyword{datasets}