savehistory()
library(oligo)

path <- "data/Affymetrix_2004-100k_trios/chip_data/Xba";
files <- list.celfiles(path=path, full.names=TRUE)[1:3]
fs <- read.celfiles(files, verbose=TRUE)
