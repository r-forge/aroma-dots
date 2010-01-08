############################################################################
# Author: Henrik Bengtsson
#
# Usage:
#  library("aroma.tcga"); library("R.menu");
#  launchMenu(aroma.tcga, "testScripts/system/");
#
# Description:
# 1. Downloads the BCRSamples.htm and stores it as a tab-delimited file.
############################################################################
if (interactive()) savehistory();
library("aroma.tcga");

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

df <- TcgaDccDownloader$readBCRSamplesTable(verbose=verbose);
str(df);

############################################################################
# HISTORY:
# 2010-01-07
# o Created.
############################################################################
