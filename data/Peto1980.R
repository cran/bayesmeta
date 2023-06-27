#
#  R. Peto. Aspirin after myocardial infarction.
#  The Lancet 315(8179):1172-1173, 1980.
#  https://doi.org/10.1016/S0140-6736(80)91626-8
#
#  P.L. Canner. An overview of six clinical trials of aspirin in coronary heart disease.
#  Statistics in Medicine, 6(3):255-263, 1987.
#  https://doi.org/10.1002/sim.4780060310
#

Peto1980 <- cbind.data.frame("publication"   =c("BrMedJ1974","JChronDis1976","Haemostasis1980",
                                                "Lancet1979","JAMA1980","Circulation1980"),
                             "study"         =c("UK-1","CDPA","GAMS","UK-2","PARIS","AMIS"),
                             "start"         =as.integer(c(1971, 1972, 1970, 1975, 1975, 1975)),
                             "end"           =as.integer(c(1973, 1975, 1977, 1979, 1979, 1979)),
                             "age"           =c(55.0, 56.5, 58.9, 56.2, 56.3, 54.8),
                             "dose"          =c(300, 972, 1500, 900, 972, 1000),
                             "followup"      =c(11.9, 22.0, 24.0, 12.0, 41.0, 39.6),
                             "treat.cases"   =as.integer(c(615, 758, 317, 832, 810, 2267)),
                             "treat.events"  =as.integer(c(49, 44, 27, 102, 85, 246)),
                             "control.cases" =as.integer(c(624, 771, 309, 850, 406, 2257)),
                             "control.events"=as.integer(c(67, 64, 32, 126, 52, 219)),
                             stringsAsFactors=FALSE)
