#
#  B.H. Rowe, J. Bretzlaff, C. Bourdon, S. Blitz, C.A. Camargo.
#  Magnesium sulfate for treating exacerbations of acute asthma
#  in the emergency department. 
#  Cochrane Database of Systematic Reviews, 1:CD001490, 2000.
#  http://dx.doi.org/10.1002/14651858.CD001490
#
#  (log-ORs for hospital admission; see page 20, Analysis 1.1)
#

RoweEtAl2000 <- cbind.data.frame("study"=c("Skobeloff (1989)", "Bloch (1995)", "Ciarallo (1997)", "Devi (1997)"),
                                 "age.range"=c("18-70 yr", "18-65 yr", "6-18 yr", "1-12 yr"),
                                 "events.MgSO4"=c(7, 7, 11, 9),
                                 "total.MgSO4"=c(19, 21, 15, 15),
                                 "events.placebo"=c(15, 11, 16, 15),
                                 "total.placebo"=c(17, 14, 16, 16),
                                 stringsAsFactors=FALSE)
