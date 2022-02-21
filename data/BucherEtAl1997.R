# Data from:
#
#   H. C. Bucher et al.
#   The results of direct and indirect treatment comparisons 
#   in meta-analysis of randomized controlled trials
#   Journal of Clinical Epidemiology, 50(6):683-691, 1997.
#   https://doi.org/10.1016/S0895-4356(97)00049-8
#
#   Table 1
#

BucherEtAl1997 <- cbind.data.frame("study"   =c("Antinori (1992)", "Mallolas (1992)", "Tocchetti (1994)",
                                                "Bozzette (1995)", "Blum (1995)", "Podzamczer (1993)",
                                                "Podzamczer (1995)", "Sirera (1995)",
                                                "Slavin (1992)", "Girard (1993)", "Torres (1993)",
                                                "Opravil (1995)", "Salmon (1995)",
                                                "Rozenbaum (1991)", "Hardy (1992)", "Schneider (1992)",
                                                "Smith (1992)", "Michelet (1993)", "May (1993)",
                                                "Stellini (1994)", "Nielsen (1995)", "Rizzardi (1995)"),
                                   "treat.A" =factor(rep(c("TMP-SMX", "AP",  "TMP-SMX"), c(8,5,9)),
                                                     levels=c("TMP-SMX", "D/P", "AP")),
                                   "treat.B" =factor(rep(c("D/P",     "D/P", "AP"),      c(8,5,9)),
                                                    levels=c("TMP-SMX", "D/P", "AP")),
                                   "events.A"=c(1, 3, 0, 42, 1, 3, 0, 6,
                                                8, 10, 15, 13, 12,
                                                0, 14, 0, 3, 1, 2, 0, 1, 5),
                                   "total.A" =c(66, 107, 15, 276, 39, 81, 104, 115,
                                                46, 176, 152, 242, 102,
                                                29, 154, 142, 27, 53, 108, 26, 47, 95),
                                   "events.B"=c(9, 8, 1, 41, 1, 13, 6, 9,
                                                9, 10, 15, 12, 5,
                                                1, 36, 6, 6, 4, 5, 2, 8, 6),
                                   "total.B" =c(63, 116, 15, 288, 47, 85, 96, 105,
                                                50, 173, 126, 291, 92,
                                                27, 156, 71, 26, 55, 106, 23, 48, 101),
                                   stringsAsFactors=FALSE)
