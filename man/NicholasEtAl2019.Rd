\name{NicholasEtAl2019}
\docType{data}
\alias{NicholasEtAl2019}
\title{Multiple sclerosis disability progression example data}
\description{Proportions of patients with disability progression in the placebo groups of 28 studies.}
\usage{data("NicholasEtAl2019")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{study}        \tab \code{character} \tab publication identifier (first author and publication year) \cr
    \bold{year}         \tab \code{numeric}   \tab publication year \cr
    \bold{patients}     \tab \code{numeric}   \tab number of placebo patients \cr
    \bold{prog.percent} \tab \code{numeric}   \tab percentage of patients with disability progression \cr
  }
}
\details{A systematic literature review investigated the characteristics of
  randomized placebo-controlled trials in multiple sclerosis published
  between 1988 and 2018 (Nicholas \emph{et al.}, 2019). A number of
  trends were observed in the trial characteristics over the
  investigated period; one of these was a decline in the proportion of
  placebo patients experiencing \emph{disability progression within 24
    months} during the course of a study. The data set contains the
  placebo groups' sizes along with the percentages of progressing
  patients within that group for 28 studies. The data were originally
  extracted from tables or Kaplan-Meier curves.
}
\source{R.S. Nicholas, E. Han, J. Raffel, J. Chataway, T. Friede.
  Over three decades study populations in progressive multiple sclerosis
  have become older and more disabled, but have lower on-trial
  progression rates: A systematic review and meta-analysis of 43
  randomised placebo-controlled trials.
  \emph{Multiple Sclerosis Journal}, \bold{25}(11):1462-1471, 2019.
  \doi{10.1177/1352458518794063}.
}
\references{C. Roever, T. Friede.
  Using the bayesmeta R package for Bayesian random-effects meta-regression.
  \emph{Computer Methods and Programs in Biomedicine},
  \bold{299}:107303, 2023.
  \doi{10.1016/j.cmpb.2022.107303}.
}
\examples{
# load data:
data("NicholasEtAl2019")

# show data:
head(NicholasEtAl2019)

\dontrun{
# compute effect sizes (logarithmic odds) from count data
# (note: effect of potential drop-outs is ignored here):
es <- escalc(measure="PLO",
             xi=patients*(prog.percent/100), ni=patients,
             slab=study, data=NicholasEtAl2019)

# illustrate estimates (log-odds):
forestplot(es, zero=NA, xlab="log(odds)", title="Nicholas et al. (2019) data")

# set up regressor matrix
# (note: "year" variable is re-scaled so that the intercept
# corresponds to the log-odds at year=2000):
X <- cbind("intercept2000" = 1, "year" = (es$year-2000))

# perform analysis:
bmr01 <- bmr(es, X=X)

# show results:
print(bmr01)
plot(bmr01)

# illustrate the data and time trend;
# first derive predictions from the model
# and specify corresponding "regressor matrix":
newx <- cbind(1, (1989:2019)-2000)

# compute credible intervals for the mean:
pred <- cbind("median"=bmr01$qpred(0.5, x=newx),
              bmr01$pred.interval(x=newx))

# compute prediction intervals:
map <- cbind("median"=bmr01$qpred(0.5, x=newx, mean=FALSE),
              bmr01$pred.interval(x=newx, mean=FALSE))

# draw empty plot:
plot(range(newx[,2]), range(map), type="n",
     xlab="publication year - 2000", ylab="log(odds)")

# show the 26 studies' estimates (and 95 percent CIs):
matlines(rbind(es$year, es$year)-2000,
         rbind(es$yi-qnorm(0.975)*sqrt(es$vi), es$yi+qnorm(0.975)*sqrt(es$vi)),
         col=1, lty=1)
points(es$year-2000, es$yi)

# show trend lines (and 95 percent intervals):
matlines(newx[,2], map, col="blue", lty=c(1,2,2))
matlines(newx[,2], pred, col="red", lty=c(1,2,2))
legend("topright", pch=15, col=c("red","blue"), c("mean","prediction"))
}
}
\keyword{datasets}
