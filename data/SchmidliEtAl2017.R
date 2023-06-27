#
#  H. Schmidli, B. Neuenschwander, T. Friede.
#  Meta-analytic-predictive use of historical variance data
#  for the design and analysis of clinical trials.
#  Computational Statistics and Data Analysis, 113:100-110, 2017.
#  https://doi.org/10.1016/j.csda.2016.08.007
#

# specify data (Table 2):
SchmidliEtAl2017 <- cbind.data.frame("study" = c("CATT", "CLEAR-IT 2", "HARBOR",
                                                 "IVAN", "VIEW 1", "VIEW 2"),
                                     "N"     = as.integer(c(599, 62, 550, 309, 909, 906)),
                                     "stdev" = c(12.11, 10.97, 10.94, 9.41, 10.97, 10.95),
                                     "df"    = as.integer(c(597, 60, 548, 307, 906, 903)),
                                     stringsAsFactors=FALSE)
