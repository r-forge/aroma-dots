library(aroma.affymetrix);
source("init.R");
source("ui.R");

dataSets <- selectDataSets(c("rawData", "probeData"));
cna(dataSets, pre=NULL);

############################################################################
# HISTORY:
# 2007-01-04
# o Created.
############################################################################
