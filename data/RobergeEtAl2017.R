# Aspirin example data
# 
#   Roberge et al.
#   The role of aspirin dose on the prevention of preeclampsia
#   and fetal growth restriction: systematic review and meta-analysis
#   American Journal of Obstetrics and Gynecology 216(2):110-120, 2017.
#   https://doi.org/10.1016/j.ajog.2016.09.076
#
# data from Table 1
# as well as Supplemental Figures 1 and 4 (preeclapsia; PE)
# and Supplemental Figures 3 and 6 (fetal growth restriction; FGR)
# FGR total number for "CLASP (1994)" study from original reference, Fig.3

RobergeEtAl2017 <- cbind.data.frame("study"     = c("Tulppala (1997)", "Benigni (1989)", "Caritis (1998a)",
                                                    "Sibai (1993a)", "Golding (1998a)", "Ebrashi (2005)",
                                                    "Zhao (2012)", "Odibo (2015)", "Porreco (1993)",
                                                    "Jamal (2012)", "Mesdaghinia (2011)", "August (1994)",
                                                    "Azar (1990)", "Bakhti (2011)", "Chiaffarino (2004)",
                                                    "Dasari (1998)", "Hermida (1997)", "Ayala (2013)",
                                                    "Michael (1992)", "Villa (2013)", "Beaufils (1985)",
                                                    "Zimmermann (1997)", "Caritis (1998b)", "CLASP (1994)",
                                                    "ECPPA (1996)", "Ferrier (1996)", "Golding (1998b)",
                                                    "Hauth (1993)", "Sibai (1993b)", "Kim (1997)",
                                                    "Wallenburg (1986)", "Wallenburg (1991)", "Byaruhanga (1998)",
                                                    "Davies (1995)", "McParland (1990)", "Rotchell (1998)",
                                                    "Wang (1996)", "Rogers (1999)", "Schrocksnadel (1992)",
                                                    "Grab (2000)", "Omrani (1992)", "Gallery (1997)",
                                                    "McCowan (1999)", "Morris (1996)", "Newnham (1995)",
                                                    "Schiff (1989)", "Trudinger (1988)", "Yu (2003)"),
                                    "year"       = c(1997, 1989, 1998, 1993, 1998, 2005, 2012, 2015, 1993,
                                                     2012, 2011, 1994, 1990, 2011, 2004, 1998, 1997, 2013,
                                                     1992, 2013, 1985, 1997, 1998, 1994, 1996, 1996, 1998,
                                                     1993, 1993, 1997, 1986, 1991, 1998, 1995, 1990, 1998,
                                                     1996, 1999, 1992, 2000, 1992, 1997, 1999, 1996, 1995,
                                                     1989, 1988, 2003),
                                    "N"          = c( 66,   33, 652,  644, 1997,  136,  237,  30,   90,
                                                      70,   80,  54,   91,   84,   35,   50, 107,  350,
                                                     110,  121,  93,   26, 1851, 2492,  606,  43, 4292,
                                                     606, 2340,  70,   46,   36,  230,  118, 100, 1485,
                                                      84,  193,  41,   43,   40,  120,   99, 102,   51,
                                                      65,   46, 554),
                                    "onset"      = factor(rep(c("up to wk 16", "after wk 16"), times=c(21,27)),
                                                          levels=c("up to wk 16", "after wk 16")),
                                    "dose"       = c(50,60,60,60,60, 75,75, 80,80,80,80, 100,100,100,100,100,100,
                                                     100,100,100, 150, 50, 60,60,60,60,60,60,60,60,60,60,
                                                     75,75,75,75,75, 80,80, 100,100,100,100,100,100,100, 150,150),
                                    "control"    = factor(c("placebo", "no treatment", "unclear")[c(1,1,1,1,1,2,1,1,1,2,2,1,2,2,2,1,1,1,1,1,2,2,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,3,1,1,1,1,1,1,1,1,1,1)],
                                                          levels=c("placebo", "no treatment", "unclear")),
                                    ############################################################################
                                    "asp.PE.events"  = c( 1,   0,   68,  22,   43,  25,   22,
                                                          3,   6,    2,   1,
                                                          3,   1,    1,  NA,   NA,   3,   11,    1,   8,  0,
                                                          4, 163,   91,  16,    1,   83,   5,   47,  NA,  0,
                                                         NA,  17,    5,   1,   10,   NA,   3,    0,
                                                          3,   2,   NA,  12,    4,   NA,   1,   NA,  49),
                                    "asp.PE.total"   = c(33,  17,  313, 320, 1009,   73, 118,
                                                         16,  48,   35,  40,
                                                         24,  46,   82,  NA,   NA,   50, 176,   55,  61, 48,
                                                         13, 941, 1259, 284,   23, 2139, 302, 1165,  NA, 23,
                                                         NA, 113,   58,  48,  739,   NA, 118,   22,
                                                         22,  21,   NA,  49,   52,   NA,  34,   NA, 276),
                                    "cont.PE.events" = c( 3,   0,   84,  24,   40,   40,  69,
                                                          3,   9,    4,   9,
                                                          5,   4,    9,  NA,   NA,    7,  22,    5,  11,  6,
                                                          2, 170,   80,  22,    1,   68,  17,   70,  NA,  7,
                                                         NA,  23,    7,  10,   12,   NA,   7,    6,
                                                          2,   7,   NA,   9,    7,   NA,   7,   NA,  52),
                                    "cont.PE.total"  = c(33,  16,  339, 324,  988,   63, 119,
                                                         14,  42,   35,  40,
                                                         25,  45,   82,  NA,   NA,   50, 174,   55,  60, 45,
                                                         13, 910, 1233, 322,   20, 2153, 302, 1175,  NA, 23,
                                                         NA, 117,   60,  52,  746,   NA,  75,   19,
                                                         21,  19,   NA,  50,   50,   NA,  31,   NA, 278),
                                    ############################################################################
                                    "asp.FGR.events" = c( 3,    2,    22,  16, NA, 13,  16,    1, NA,  1,  0,
                                                          0,   NA,     0,   2,  1,  1,  16,   NA,  2,  4,
                                                          2,  111,    93,  28, NA, NA,  17,   61,  5,  4,
                                                         NA,   18,     6,   7, NA,  3,
                                                         NA,    1,    NA,  NA, NA, 37,  14,   25,  2, NA,  61),
                                    "asp.FGR.total"  = c(23,   17,   329, 318, NA, 73, 118,   16, NA, 35, 40,
                                                         24,   NA,    82,  16, 25, 50, 176,   NA, 61, 48,
                                                         13, 1255,  1321, 286, NA, NA, 302, 1165, 32, 23,
                                                         NA,  114,    58,  48, NA, 40,
                                                         NA,   22,    NA,  NA, NA, 49,  52,   29, 34, NA, 276),
                                    "cont.FGR.events"= c( 3,    6,    28,  21, NA, 21,  36,    1, NA,  2,  0,
                                                          1,   NA,     1,   5,  5,  2,  32,   NA,  6, 13,
                                                          1,   85,    87,  42, NA, NA,  19,   76,  5,  6,
                                                         NA,   20,     6,   7, NA, 12,
                                                         NA,    2,    NA,  NA, NA, 39,  11,   27,  6, NA,  68),
                                    "cont.FGR.total" = c(23,   16,   374, 324, NA, 63, 119,   14, NA, 35, 40,
                                                         25,   NA,    82,  19, 25, 50, 174,   NA, 60, 45,
                                                         13, 1182,  1301, 329, NA, NA, 302, 1180, 38, 23,
                                                         NA,  122,    60,  52, NA, 44,
                                                         NA,   19,    NA,  NA, NA, 50,  50,   30, 31, NA, 278),
                                    stringsAsFactors=FALSE)

# check for missings:
nodata <- which(!(complete.cases(RobergeEtAl2017[,7:10]) | (complete.cases(RobergeEtAl2017[,11:14]))))

# drop missings:
RobergeEtAl2017 <- RobergeEtAl2017[-nodata,]
rm(list="nodata")

# re-assign row names (-numbers):
rownames(RobergeEtAl2017) <- sprintf("%02d", 1:nrow(RobergeEtAl2017))
