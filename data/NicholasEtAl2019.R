#
# EDSS progression data from
#
#   R. S. Nicholas, E. Han, J. Raffel, J. Chataway, T. Friede.
#   Over three decades study populations in progressive multiple sclerosis
#   have become older and more disabled, but have lower on-trial progression rates:
#   A systematic review and meta-analysis of 43 randomised placebo-controlled trials.
#   Multiple Sclerosis Journal, 25(11):1462-1471, 2019.
#   https://doi.org/10.1177/1352458518794063
#
# (Fig. 7).
#

NicholasEtAl2019 <- cbind.data.frame("study"       =c("Kastrukoff (1990)", "Wolinsky (1990)", "Bornstein (1991)",
                                                      "Likosky (1991)", "Noseworthy (1991)", "Milanese (1992)",
                                                      "Milligan (1994)", "Goodkin (1995)", "Kappos (1998)",
                                                      "Francis (2001)", "Cohen (2002)", "Hartung (2002)",
                                                      "Leary (2003)", "Andersen (2004)", "Hommes (2004)",
                                                      "Panitch (2004)", "Warren (2006)", "Pohlau (2007)",
                                                      "Wolinsky (2007)", "Hawker (2009)", "Montalban (2009)",
                                                      "Freedman (2011)", "Mostert (2013)", "Zajicek (2013)",
                                                      "Chataway (2014)", "Lublin (2016)", "Montalban (2016)",
                                                      "Kapoor (2018)"),
                                     "year"        =c(1990, 1990, 1991,  1991, 1991, 1992,  1994, 1995, 1998,
                                                      2001, 2002, 2002,  2003, 2004, 2004,  2004, 2006, 2007,
                                                      2007, 2009, 2009,  2011, 2013, 2013,  2014, 2016, 2016,
                                                      2018),
                                     "patients"    =c( 50, 274, 55,   20,  56,  21,   18,  29, 358,
                                                      205, 219, 64,   20, 178, 159,  308,  16, 115,
                                                      316, 147, 37,  305,  22, 164,   70, 487, 244,
                                                      449),
                                     "prog.percent"=c(46.0, 72.5, 29.5,  80.0, 22.0, 55.0,  61.1, 51.7, 48.0,
                                                      53.0, 33.7, 22.0,  43.0, 34.0, 44.0,  28.0, 56.0, 63.0,
                                                      32.0, 38.5, 43.2,  29.2, 32.0, 28.0,  54.0, 43.0, 31.0,
                                                      15.0),
                                     stringsAsFactors=FALSE)
