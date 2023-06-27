#
#   B. Neuenschwander, S. Weber, H. Schmidli, A. O'Hagan.
#   Predictively consistent prior effective sample sizes.
#   Biometrics 76(2):578-587, 2020.
#   https://doi.org/10.1111/biom.13252
#
#   D. Baeten et al.
#   Anti-interleukin-17A monoclonal antibody secukinumab
#   in treatment of ankylosing spondylitis: a randomised,
#   double-blind, placebo-controlled trial.
#   The Lancet 382(9906):1705-1713, 2013.
#   https://doi.org/10.1016/S0140-6736(13)61134-4
#

BaetenEtAl2013 <- cbind.data.frame("study"  = c("ATLAS", "Canadian AS", "Wyeth", "Calin",
                                                "Davis", "Gorman", "ASSERT", "Braun"),
                                   "year"   = c(2005, 2005, 2006, 2003,
                                                2003, 2002, 2005, 2002),
                                   "events" = c( 23, 12, 19,  9,   39,  6,  9, 10),
                                   "total"  = c(107, 44, 51, 39,  139, 20, 78, 35),
                                   stringsAsFactors=FALSE)

rownames(BaetenEtAl2013) <- paste0("study_",1:8)
