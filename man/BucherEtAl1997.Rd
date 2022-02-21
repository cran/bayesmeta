\name{BucherEtAl1997}
\docType{data}
\alias{BucherEtAl1997}
\title{Direct and indirect comparison example data}
\description{Numbers of subjects and events in the different treatment arms of 22 studies.}
\usage{data("BucherEtAl1997")}
\format{The data frame contains the following columns:
  \tabular{lll}{
    \bold{study}        \tab \code{character} \tab publication identifier (first author and publication year) \cr
    \bold{treat.A}      \tab \code{factor}    \tab treatment in first study arm (\dQuote{TMP-SMX} or \dQuote{AP}) \cr
    \bold{treat.B}      \tab \code{factor}    \tab treatment in second study arm (\dQuote{D/P} or \dQuote{AP}) \cr
    \bold{events.A}     \tab \code{numeric}   \tab number of events in first study arm \cr
    \bold{events.B}     \tab \code{numeric}   \tab number of events in second study arm \cr
    \bold{total.A}      \tab \code{numeric}   \tab total number of patients in first study arm \cr
    \bold{total.B}      \tab \code{numeric}   \tab total number of patients in second study arm \cr
  }
}
\details{Bucher \emph{et al.} (1997) discussed the example case of the
  comparison of \emph{sulphametoxazole-trimethoprim (TMP-SMX)} versus
  \emph{dapsone/pyrimethamine (D/P)} for the prophylaxis of
  \emph{Pneumocystis carinii} pneumonia in HIV patients. Eight
  studies had undertaken a head-to-head comparison of both medications,
  but an additional 14 studies were available investigating one of the
  two medications with \emph{aerosolized pentamidine (AP)} as a
  comparator. Nine studies compared TMP-SMX vs. AP, and five studies
  compared D/P vs. AP. Together these provide \emph{indirect} evidence
  on the effect of TMP-SMX compared to D/P (Kiefer \emph{et al.}, 2015).

  The example constitutes a simple case of a \emph{network meta-analysis
  (NMA)} setup, where only two-armed studies are considered, and
  analysis is based on pairwise comparisons of treatments (or
  \emph{contrasts}). In this case, the joint analysis of \emph{direct}
  and \emph{indirect} evidence may be implemented as a special case of a
  meta-regression (Higgins \emph{et al.}, 2019; Sec. 11.4.2).
  The original data in fact included some three-armed studies, in which
  case one of the arms was deliberately omitted (Bucher \emph{et al.}; 1997).
}
\source{H.C. Bucher, G.H. Guyatt, L.E. Griffith, S.D. Walter.
  The results of direct and indirect treatment comparisons 
  in meta-analysis of randomized controlled trials.
  \emph{Journal of Clinical Epidemiology}, \bold{50}(6):683-691, 1997.
  \doi{10.1016/S0895-4356(97)00049-8}.
}
\references{J.P.T. Higgins, J. Thomas, J. Chandler, M. Cumpston, T. Li,
  M.J. Page, V.A. Welch (eds.).
  \emph{Cochrane handbook for systematic reviews of interventions}.
  Wiley and Sons, 2nd edition, 2019.
  \doi{10.1002/9781119536604}.
  \url{http://training.cochrane.org/handbook}.

  C. Kiefer, S. Sturtz, R. Bender.
  Indirect comparisons and network meta-analyses.
  \emph{Deutsches Aerzteblatt International},
  \bold{112}(47):803-808, 2015.
  \doi{10.3238/arztebl.2015.0803}.
}
\examples{
# load data:
data("BucherEtAl1997")

# show data:
head(BucherEtAl1997)

\dontrun{
# compute effect sizes (log-ORs for pairwise comparisons)
# from the count data:
es <- escalc(measure="OR",
             ai=events.A, n1i=total.A,   # "exposure group"
             ci=events.B, n2i=total.B,   # "control group"
             slab=study, data=BucherEtAl1997)

# specify regressor matrix:
X <- cbind("TMP.DP" = rep(c(1, 0, 1), c(8,5,9)),
           "AP.DP"  = rep(c(0, 1,-1), c(8,5,9)))

# perform Bayesian meta-regression:
bmr01 <- bmr(es, X=X)

# show default output:
print(bmr01)

# specify contrast matrix:
contrastX <- rbind("TMP-SMX vs. D/P"=c(1,0),
                   "AP vs. D/P"     =c(0,1),
                   "TMP-SMX vs. AP" =c(1,-1))
# show summary including contrast estimates:
summary(bmr01, X.mean=contrastX)
# show forest plot including contrast estimates:
forestplot(bmr01, X.mean=contrastX, xlab="log-OR")


# perform frequentist meta-regression:
fmr01 <- rma(es, mods=X, intercept=FALSE)
print(fmr01)

# compare Bayesian and frequentist results;
# estimated log-OR for "TMP-SMX" vs. "D/P"
rbind("bayesmeta"=bmr01$summary[c("mean","sd"),"TMP.DP"],
      "rma"      =c(fmr01$beta["TMP.DP",], fmr01$se[1]))

# estimated log-OR for "AP" vs. "D/P"
rbind("bayesmeta"=bmr01$summary[c("mean","sd"),"AP.DP"],
      "rma"      =c(fmr01$beta["AP.DP",], fmr01$se[2]))

# estimated heterogeneity:
rbind("bayesmeta"=bmr01$summary["median","tau"],
      "rma"      =sqrt(fmr01$tau2))
}
}
\keyword{datasets}
