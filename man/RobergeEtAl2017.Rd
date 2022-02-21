\name{RobergeEtAl2017}
\docType{data}
\alias{RobergeEtAl2017}
\title{Aspirin during pregnancy example data}
\description{Numbers of cases (patients) and events (preeclampsia (PE) or
  fetal growth restriction (FGR)) in experimental and control groups of
  45 studies.}
\usage{data("RobergeEtAl2017")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{study}           \tab \code{character} \tab publication identifier (first author and publication year) \cr
    \bold{year}            \tab \code{numeric}   \tab publication year \cr
    \bold{N}               \tab \code{numeric}   \tab number of patients \cr
    \bold{onset}           \tab \code{factor}    \tab treatment onset (up to 16 weeks' gestation or later) \cr
    \bold{dose}            \tab \code{numeric}   \tab dose (mg/day) \cr
    \bold{control}         \tab \code{factor}    \tab type of control group \cr
    \bold{asp.PE.events}   \tab \code{numeric}   \tab number of PE events in aspirin group \cr
    \bold{asp.PE.total}    \tab \code{numeric}   \tab number of PE cases in aspirin group \cr
    \bold{cont.PE.events}  \tab \code{numeric}   \tab number of PE events in control group \cr
    \bold{cont.PE.total}   \tab \code{numeric}   \tab number of PE cases in control group \cr
    \bold{asp.FGR.events}  \tab \code{numeric}   \tab number of FGR events in aspirin group \cr
    \bold{asp.FGR.total}   \tab \code{numeric}   \tab number of FGR cases in aspirin group \cr
    \bold{cont.FGR.events} \tab \code{numeric}   \tab number of FGR events in control group \cr
    \bold{cont.FGR.total}  \tab \code{numeric}   \tab number of FGR cases in control group \cr
  }
}
\details{A systematic literature review was performed in order to
  summarize the evidence on effects of aspirin administered during
  pregnancy. Of particular interest were occurrences of
  \emph{preeclampsia (PE)} and \emph{fetal growth restriction (FGR)}. A
  total of 45 relevant randomized controlled trials (RCTs) were found,
  out of which 40 reported on PE, and 35 reported on FGR. Besides event
  rates, the mode of administration (treatment onset (early vs. late)
  and dose (in mg)) was also recorded for each study.
}
\source{S. Roberge, K. Nicolaides, S. Demers, J. Hyett, N. Chaillet, E. Bujold.
  The role of aspirin dose on the prevention of preeclampsia and fetal
  growth restriction: systematic review and meta-analysis.
  \emph{American Journal of Obstetrics & Gynecology},
  \bold{216}(2):110-120, 2017.
  \doi{10.1016/j.ajog.2016.09.076}.
}
\seealso{
  \code{\link{bmr}}, \code{\link[metafor]{escalc}},
  \code{\link[stats]{model.matrix}}.
}
\examples{
# load data:
data("RobergeEtAl2017")
str(RobergeEtAl2017)
head(RobergeEtAl2017)

# compute effect sizes (log odds ratios) from count data
# (using the "metafor" package's "escalc()" function);
# preeclampsia (PE):
es.pe  <- escalc(measure="OR",
                 ai=asp.PE.events,  n1i=asp.PE.total,
                 ci=cont.PE.events, n2i=cont.PE.total,
                 slab=study, data=RobergeEtAl2017,
                 subset=complete.cases(RobergeEtAl2017[,7:10]))
# show forest plot:
forestplot(es.pe, title="preeclampsia (PE)")
# show "bubble plot" (bubble sizes are
# inversely proportional to standard errors):
plot(es.pe$dose, es.pe$yi, cex=1/sqrt(es.pe$vi),
     col=c("blue","red")[as.numeric(es.pe$onset)],
     xlab="dose (mg)", ylab="log-OR (PE)", main="Roberge et al. (2017)")
legend("topright", col=c("blue","red"), c("early onset", "late onset"), pch=1)

# fetal growth restriction (FGR):
es.fgr <- escalc(measure="OR",
                 ai=asp.FGR.events,  n1i=asp.FGR.total,
                 ci=cont.FGR.events, n2i=cont.FGR.total,
                 slab=study, data=RobergeEtAl2017,
                 subset=complete.cases(RobergeEtAl2017[,11:14]))
# show forest plot:
forestplot(es.fgr, title="fetal growth restriction (FGR)")
# show "bubble plot":
plot(es.fgr$dose, es.fgr$yi, cex=1/sqrt(es.fgr$vi),
     col=c("blue","red")[as.numeric(es.fgr$onset)],
     xlab="dose (mg)", ylab="log-OR (FGR)", main="Roberge et al. (2017)")
legend("topright", col=c("blue","red"), c("early onset", "late onset"), pch=1)

\dontrun{
# set up regressor matrix (common intercept and slope):
X01 <- model.matrix(~ dose, data=es.pe)
colnames(X01) <- c("intercept", "slope")
print(X01)

# perform regression:
bmr01 <- bmr(es.pe, X=X01)
bmr01$summary

# set up alternative regressor matrix
# (individual intercepts and slopes for two subgroups):
X02 <- model.matrix(~ -1 + onset + onset:dose, data=es.pe)
colnames(X02) <- c("intEarly", "intLate", "slopeEarly", "slopeLate")
print(X02)

# perform regression:
bmr02 <- bmr(es.pe, X=X02)
bmr02$summary

# derive predictions from the model;
# specify corresponding "regressor matrices":
newx.early <- cbind(1, 0, seq(50, 150, by=5), 0)
newx.late  <- cbind(0, 1, 0, seq(50, 150, by=5))
# (note: columns correspond to "beta" parameters)

# compute predicted medians and 95 percent intervals: 
pred.early <- cbind("median"=bmr02$qpred(0.5, x=newx.early),
                    bmr02$pred.interval(x=newx.early))
pred.late <- cbind("median"=bmr02$qpred(0.5, x=newx.late),
                    bmr02$pred.interval(x=newx.late))

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
\keyword{datasets}
